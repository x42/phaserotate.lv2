/*
 * Copyright (C) 2021 Robin Gareus <robin@gareus.org>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#define PHASEROTATE_URI "http://gareus.org/oss/lv2/phaserotate"

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include <lv2/lv2plug.in/ns/ext/log/logger.h>
#include <lv2/lv2plug.in/ns/lv2core/lv2.h>

static pthread_mutex_t fftw_planner_lock = PTHREAD_MUTEX_INITIALIZER;
static unsigned int    instance_count    = 0;

typedef struct {
	/* ports */
	float* p_in;
	float* p_out;
	float* p_angle;
	float* p_latency;

	/* Config */
	uint32_t fftlen;
	uint32_t parsiz; // == fftlen / 2

	/* Latent Buffers */
	float*   buf_src; // input to be processed by FFT
	float*   buf_out; // result of inv FFT
	uint32_t offset;  // r/w offset in latent buffers
	uint32_t overlap; // input overlap offset (A/B; 0 or parsiz)

	/* FFT Bufers and Plan */
	float*         window;
	float*         time_data;
	fftwf_complex* freq_data;
	fftwf_plan     plan_r2c;
	fftwf_plan     plan_c2r;

} FFTiProc;

/* *****************************************************************************
 * LV2 Plugin
 */

static void
cleanup (LV2_Handle instance)
{
	FFTiProc* self = (FFTiProc*)instance;
	free (self->window);
	free (self->buf_src);
	free (self->buf_out);
	fftwf_free (self->time_data);
	fftwf_free (self->freq_data);

	pthread_mutex_lock (&fftw_planner_lock);
	fftwf_destroy_plan (self->plan_r2c);
	fftwf_destroy_plan (self->plan_c2r);

	if (instance_count > 0) {
		--instance_count;
	}

#ifdef WITH_STATIC_FFTW_CLEANUP
	/* use this only when statically linking to a local fftw!
	 *
	 * "After calling fftw_cleanup, all existing plans become undefined,
	 *  and you should not attempt to execute them nor to destroy them."
	 * [http://www.fftw.org/fftw3_doc/Using-Plans.html]
	 *
	 * If libfftwf is shared with other plugins or the host this can
	 * cause undefined behavior.
	 */
	if (instance_count == 0) {
		fftwf_cleanup ();
	}
#endif
	pthread_mutex_unlock (&fftw_planner_lock);

	free (instance);
}

static LV2_Handle
instantiate (const LV2_Descriptor*     descriptor,
             double                    rate,
             const char*               bundle_path,
             const LV2_Feature* const* features)
{
	FFTiProc* self = (FFTiProc*)calloc (1, sizeof (FFTiProc));

	if (rate < 64000) {
		/* 44.1 and 48 kHz */
		self->fftlen = 2048;
	} else if (rate < 128000) {
		/* 88.2 or 96 kHz. */
		self->fftlen = 4096;
	} else {
		self->fftlen = 8192;
	}

	self->parsiz    = self->fftlen / 2;
	self->buf_src   = (float*)calloc (self->fftlen, sizeof (float));
	self->buf_out   = (float*)calloc (self->parsiz, sizeof (float));
	self->window    = (float*)calloc (self->fftlen, sizeof (float));
	self->time_data = (float*)fftwf_malloc (self->fftlen * sizeof (float));
	self->freq_data = (fftwf_complex*)fftwf_malloc ((self->parsiz + 1) * sizeof (fftwf_complex));

	if (!self->buf_src || !self->buf_out || !self->window || !self->time_data || !self->freq_data) {
		pthread_mutex_lock (&fftw_planner_lock);
		++instance_count;
		pthread_mutex_unlock (&fftw_planner_lock);
		cleanup ((LV2_Handle)self);
		return NULL;
	}

	pthread_mutex_lock (&fftw_planner_lock);
	++instance_count;
	self->plan_r2c = fftwf_plan_dft_r2c_1d (self->fftlen, self->time_data, self->freq_data, FFTW_ESTIMATE);
	self->plan_c2r = fftwf_plan_dft_c2r_1d (self->fftlen, self->freq_data, self->time_data, FFTW_ESTIMATE);
	pthread_mutex_unlock (&fftw_planner_lock);

	if (!self->plan_r2c || !self->plan_c2r) {
		cleanup ((LV2_Handle)self);
		return NULL;
	}

	/* create normalized raised cosine window */
	const double norm = 0.25 / self->parsiz;
	for (uint32_t i = 0; i < self->fftlen; ++i) {
		self->window[i] = norm * (1.0 - cos (2.0 * M_PI * i / self->fftlen));
	}

	self->overlap = 0;

	return (LV2_Handle)self;
}

static void
connect_port (LV2_Handle instance,
              uint32_t   port,
              void*      data)
{
	FFTiProc* self = (FFTiProc*)instance;

	switch (port) {
		case 0:
			self->p_in = (float*)data;
			break;
		case 1:
			self->p_out = (float*)data;
			break;
		case 2:
			self->p_angle = (float*)data;
			break;
		case 3:
			self->p_latency = (float*)data;
			break;
		default:
			break;
	}
}

static void
activate (LV2_Handle instance)
{
	FFTiProc* self = (FFTiProc*)instance;
	memset (self->buf_src, 0, self->fftlen * sizeof (float));
	memset (self->buf_out, 0, self->parsiz * sizeof (float));
	memset (self->time_data, 0, self->fftlen * sizeof (float));
}

static void
run (LV2_Handle instance, uint32_t n_samples)
{
	FFTiProc* self  = (FFTiProc*)instance;
	double    angle = *self->p_angle / 360.0;

	if (angle < 0) {
		angle = 0;
	}
	if (angle > 1.0) {
		angle = 1.0;
	}

	const float ca = cos (2.0 * M_PI * angle);
	const float sa = sin (2.0 * M_PI * angle);

	/* copy/forward no-inplace buffers */
	if (self->p_in != self->p_out) {
		memcpy (self->p_out, self->p_in, sizeof (float) * n_samples);
	}

	/* process output buffer in-place */
	float* iobuf = self->p_out;

	/* localize constants */
	const uint32_t parsiz = self->parsiz;
	const uint32_t fftlen = self->fftlen;

	/* localize state */
	uint32_t remain  = n_samples;
	uint32_t offset  = self->offset;
	uint32_t overlap = self->overlap;

	while (remain > 0) {
		uint32_t ns = parsiz - offset;
		if (ns > remain) {
			ns = remain;
		}

		/* copy latent buffers */
		memcpy (&self->buf_src[offset + overlap], iobuf, sizeof (float) * ns);
		memcpy (iobuf, &self->buf_out[offset], sizeof (float) * ns);

		iobuf  += ns;
		offset += ns;
		remain -= ns;

		if (offset == parsiz) {
			offset = 0;

			/* copy end/overlap of prev. iFFT */
			memcpy (self->buf_out, self->time_data + parsiz, sizeof (float) * parsiz);

			/* fill fft buffer */
			if (0 != overlap) {
				memcpy (self->time_data, self->buf_src, sizeof (float) * parsiz);
				memcpy (self->time_data + parsiz, self->buf_src + parsiz, sizeof (float) * parsiz);
				overlap = 0;
			} else {
				memcpy (self->time_data, self->buf_src + parsiz, sizeof (float) * parsiz);
				memcpy (self->time_data + parsiz, self->buf_src, sizeof (float) * parsiz);
				overlap = parsiz;
			}

#pragma GCC ivdep /* do not check for aliasing, buffers do not overlap */
			/* apply window */
			for (uint32_t k = 0; k < fftlen; ++k) {
				self->time_data[k] *= self->window[k];
			}

			fftwf_execute_dft_r2c (self->plan_r2c, self->time_data, self->freq_data);

			/* rotate phase */
			for (uint32_t k = 0; k <= parsiz; ++k) {
				const float re        = self->freq_data[k][0];
				self->freq_data[k][0] = (ca * re - sa * self->freq_data[k][1]);
				self->freq_data[k][1] = (sa * re + ca * self->freq_data[k][1]);
			}

			fftwf_execute_dft_c2r (self->plan_c2r, self->freq_data, self->time_data);

#pragma GCC ivdep
			for (uint32_t k = 0; k < parsiz; ++k) {
				self->buf_out[k] += self->time_data[k];
			}
		}
	}

	/* copy back state */
	self->overlap = overlap;
	self->offset  = offset;

	/* announce latency */
	*self->p_latency = self->fftlen;
}

static const void*
extension_data (const char* uri)
{
	return NULL;
}

static const LV2_Descriptor descriptor = {
	PHASEROTATE_URI,
	instantiate,
	connect_port,
	activate,
	run,
	NULL,
	cleanup,
	extension_data
};

/* clang-format off */
#undef LV2_SYMBOL_EXPORT
#ifdef _WIN32
# define LV2_SYMBOL_EXPORT __declspec(dllexport)
#else
# define LV2_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
#endif
/* clang-format on */
LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor (uint32_t index)
{
	switch (index) {
		case 0:
			return &descriptor;
		default:
			return NULL;
	}
}
