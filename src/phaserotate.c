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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE // needed for M_PI
#endif

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

	/* Parameter */
	float angle;
	float sa; // sin (angle)
	float ca; // cos (angle)

	/* Config */
	uint32_t fftlen;
	uint32_t firlen;
	uint32_t parsiz; // fftlen / 2
	uint32_t firlat; // firlen / 2
	uint32_t n_segm; // firlen / fftlen

	float interp_th; // parsiz / 1e6
	float interp_nm; // 1 / parsiz

	/* Latent Buffers */
	float*   buf_src; // input to be processed by FFT
	float*   buf_out; // result of inv FFT
	uint32_t offset;  // r/w offset in latent buffers
	int32_t  overlap; // input overlap offset

	/* FFT Bufers and Plan */
	float*          time_data;
	fftwf_complex*  freq_data;
	fftwf_complex*  freq_sum;
	fftwf_complex** freq_fir;
	fftwf_plan      plan_r2c;
	fftwf_plan      plan_c2r;

} FFTiProc;

static void
sin_cos (float angle, float* s, float* c)
{
	static const float twopi = 2 * M_PI;
#ifdef __APPLE__
	angle *= twopi;
	*s = sinf (angle);
	*c = cosf (angle);
#else
	sincosf (angle * twopi, s, c);
#endif
}

static void
cleanup (LV2_Handle instance)
{
	FFTiProc* self = (FFTiProc*)instance;

	fftwf_free (self->time_data);
	fftwf_free (self->freq_data);
	fftwf_free (self->freq_sum);

	if (self->freq_fir) {
		for (uint32_t i = 0; i < self->n_segm; ++i) {
			fftwf_free (self->freq_fir[i]);
		}
	}
	free (self->freq_fir);
	free (self->buf_src);
	free (self->buf_out);

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
		self->fftlen = 512;
		self->firlen = 3072;
	} else if (rate < 128000) {
		/* 88.2 or 96 kHz. */
		self->fftlen = 1024;
		self->firlen = 4096;
	} else {
		self->fftlen = 2048;
		self->firlen = 8192;
	}

	self->angle     = 0;
	self->overlap   = 0;
	self->parsiz    = self->fftlen / 2;
	self->firlat    = self->firlen / 2;
	self->n_segm    = self->firlen / self->parsiz;
	self->interp_th = self->parsiz * 1e-6f;
	self->interp_nm = 1.f / self->parsiz;

	self->buf_src   = (float*)calloc (self->firlen, sizeof (float));
	self->buf_out   = (float*)calloc (self->parsiz, sizeof (float));
	self->time_data = (float*)fftwf_malloc (self->fftlen * sizeof (float));
	self->freq_data = (fftwf_complex*)fftwf_malloc ((self->parsiz + 1) * sizeof (fftwf_complex));
	self->freq_sum  = (fftwf_complex*)fftwf_malloc ((self->parsiz + 1) * sizeof (fftwf_complex));

	if (!self->buf_src || !self->buf_out || !self->time_data || !self->freq_data || !self->freq_sum) {
		pthread_mutex_lock (&fftw_planner_lock);
		++instance_count;
		pthread_mutex_unlock (&fftw_planner_lock);
		cleanup ((LV2_Handle)self);
		return NULL;
	}

	self->freq_fir = (fftwf_complex**)calloc (self->n_segm, sizeof (fftwf_complex*));

	if (!self->freq_fir) {
		pthread_mutex_lock (&fftw_planner_lock);
		++instance_count;
		pthread_mutex_unlock (&fftw_planner_lock);
		cleanup ((LV2_Handle)self);
		return NULL;
	}
	for (uint32_t i = 0; i < self->n_segm; ++i) {
		self->freq_fir[i] = (fftwf_complex*)fftwf_malloc ((self->parsiz + 1) * sizeof (fftwf_complex));
	}
	for (uint32_t i = 0; i < self->n_segm; ++i) {
		if (!self->freq_fir[i]) {
			pthread_mutex_lock (&fftw_planner_lock);
			++instance_count;
			pthread_mutex_unlock (&fftw_planner_lock);
			cleanup ((LV2_Handle)self);
			return NULL;
		}
	}

	float*         fir          = (float*)fftwf_malloc (self->firlen * sizeof (float));
	fftwf_complex* freq_hilbert = (fftwf_complex*)fftwf_malloc ((self->firlat + 1) * sizeof (fftwf_complex));

	if (!fir || !freq_hilbert) {
		fftwf_free (fir);
		fftwf_free (freq_hilbert);
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

	fftwf_plan plan_fir_c2r = fftwf_plan_dft_c2r_1d (self->firlen, freq_hilbert, fir, FFTW_ESTIMATE);
	pthread_mutex_unlock (&fftw_planner_lock);

	if (!self->plan_r2c || !self->plan_c2r || !plan_fir_c2r) {
		fftwf_free (fir);
		fftwf_free (freq_hilbert);
		cleanup ((LV2_Handle)self);
		return NULL;
	}

	/* generate hilbert FIR */
	for (uint32_t i = 0; i <= self->firlat; ++i) {
		const float im     = i & 1 ? -1.f : 1.f;
		freq_hilbert[i][0] = 0;
		freq_hilbert[i][1] = im;
	}
	fftwf_execute_dft_c2r (plan_fir_c2r, freq_hilbert, fir);

	pthread_mutex_lock (&fftw_planner_lock);
	fftwf_destroy_plan (plan_fir_c2r);
	pthread_mutex_unlock (&fftw_planner_lock);

	/* normalize and segment FIR */
	const double fnorm = 0.5 / self->firlen;
	const double flen_ = 1.0 / self->firlen;
	for (uint32_t i = 0; i < self->firlen; ++i) {
		fir[i] *= fnorm * (1 - cos (2.0 * M_PI * i * flen_));
	}

	const float norm = 0.5 / self->parsiz;
	memset (&self->time_data[self->parsiz], 0, sizeof (float) * self->parsiz);

	for (uint32_t s = 0; s < self->n_segm; ++s) {
		for (int i = 0; i < self->parsiz; ++i) {
			self->time_data[i] = norm * fir[s * self->parsiz + i];
		}
		fftwf_execute_dft_r2c (self->plan_r2c, self->time_data, self->freq_fir[s]);
	}

	fftwf_free (freq_hilbert);
	fftwf_free (fir);

	sin_cos (self->angle, &self->sa, &self->ca);

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
	memset (self->buf_src, 0, self->firlen * sizeof (float));
	memset (self->buf_out, 0, self->parsiz * sizeof (float));
	memset (self->time_data, 0, self->fftlen * sizeof (float));
	self->overlap = 0;
	self->offset  = 0;
}

static void
run (LV2_Handle instance, uint32_t n_samples)
{
	FFTiProc* self = (FFTiProc*)instance;

	float angle        = self->angle;
	float target_angle = *self->p_angle / -360.f;

	if (target_angle < -.5f) {
		target_angle = -.5f;
	}
	if (target_angle > 0.5f) {
		target_angle = 0.5f;
	}

	/* copy/forward no-inplace buffers */
	if (self->p_in != self->p_out) {
		memcpy (self->p_out, self->p_in, sizeof (float) * n_samples);
	}

	/* process output buffer in-place */
	float* iobuf = self->p_out;

	/* localize constants */
	const uint32_t parsiz = self->parsiz;
	const uint32_t firlen = self->firlen;
	const uint32_t firlat = self->firlat;
	const uint32_t n_segm = self->n_segm;

	/* localize FFT buffer pointers */
	float*                      time_data = self->time_data;
	fftwf_complex*              freq_data = self->freq_data;
	fftwf_complex*              freq_sum  = self->freq_sum;
	fftwf_complex* const* const freq_fir  = self->freq_fir;

	/* localize state */
	uint32_t remain  = n_samples;
	uint32_t offset  = self->offset;
	int32_t  overlap = self->overlap;

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
			memcpy (self->buf_out, &time_data[parsiz], sizeof (float) * parsiz);

			/* clear convolution buffers */
			memset (&time_data[parsiz], 0, sizeof (float) * parsiz);
			memset (freq_sum, 0, sizeof (fftwf_complex) * (parsiz + 1));

			int32_t olp = overlap;
			for (uint32_t s = 0; s < n_segm; ++s) {
				memcpy (time_data, &self->buf_src[olp], sizeof (float) * parsiz);
				olp -= parsiz;
				if (olp < 0) {
					olp += firlen;
				}

				fftwf_complex* ffir = freq_fir[s];
				fftwf_execute_dft_r2c (self->plan_r2c, time_data, freq_data);
				/* FIR convolution, 90deg phase shift */
#pragma GCC ivdep /* do not check for aliasing, buffers do not overlap */
				for (uint32_t i = 0; i <= parsiz; ++i) {
					freq_sum[i][0] += freq_data[i][0] * ffir[i][0] - freq_data[i][1] * ffir[i][1];
					freq_sum[i][1] += freq_data[i][0] * ffir[i][1] + freq_data[i][1] * ffir[i][0];
				}
			}

			fftwf_execute_dft_c2r (self->plan_c2r, freq_sum, time_data);

#pragma GCC ivdep
			for (uint32_t i = 0; i < parsiz; ++i) {
				self->buf_out[i] += time_data[i];
			}

			/* delayed input (firlen / 2) */
			int32_t off = overlap - firlat;
			if (off < 0) {
				off += firlen;
			}

			float* in  = &self->buf_src[off];
			float* out = self->buf_out;

			if (target_angle != angle) {
				/* interpolate */
				float da = (target_angle - angle);
				if (fabs (da) > 0.5) {
					/* wrap around at +/- 180 */
					if (da < 0) {
						da += 1.f;
					} else {
						da -= 1.f;
					}
				}

				da *= self->interp_nm;

				int         final  = 0;
				const float thresh = self->interp_th;
				if (da > thresh) {
					da = thresh;
				} else if (da < -thresh) {
					da = -thresh;
				} else {
					final = 1;
				}

				for (uint32_t i = 0; i < parsiz; ++i) {
					float ca, sa;
					sin_cos (angle, &sa, &ca);
					out[i] = ca * in[i] + sa * out[i];
					angle += da;
				}

				if (final) {
					angle = target_angle;
				}
				if (angle == target_angle) {
					sin_cos (angle, &self->sa, &self->ca);
				}
			} else {
				const float sa = self->sa;
				const float ca = self->ca;
#pragma GCC ivdep
				for (uint32_t i = 0; i < parsiz; ++i) {
					out[i] = ca * in[i] + sa * out[i];
				}
			}

			overlap = (overlap + parsiz) % firlen;
		}
	}

	/* copy back state */
	self->overlap = overlap;
	self->offset  = offset;
	self->angle   = angle;

	/* announce latency */
	*self->p_latency = parsiz + firlat;
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
