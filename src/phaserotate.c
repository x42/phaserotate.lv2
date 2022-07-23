/*
 * Copyright (C) 2021,2022 Robin Gareus <robin@gareus.org>
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

#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "lv2/lv2plug.in/ns/ext/options/options.h"
#include <lv2/lv2plug.in/ns/ext/log/logger.h>
#include <lv2/lv2plug.in/ns/lv2core/lv2.h>

#include "phaserotate.h"

static pthread_mutex_t fftw_planner_lock = PTHREAD_MUTEX_INITIALIZER;
static unsigned int    instance_count    = 0;

typedef struct {
	/* Ports */
	float* p_in;
	float* p_out;
	float* p_angle;

	/* Parameter */
	float angle;
	float sa; // sin (angle)
	float ca; // cos (angle)

	/* Latent Buffers */
	float*   buf_dly; // delayed input for input meter
	float*   buf_src; // input to be processed by FFT
	float*   buf_out; // result of inv FFT
	uint32_t offset;  // r/w offset in latent buffers
	int32_t  overlap; // input overlap offset

	/* FFT Buffer*/
	float* time_data;

	/* Meter */
	int   mtr_reset_delay;
	int   mtr_holdcnt[2];
	float mtr_momentary[2];
	float mtr_peak[2];
	float mtr_diff[2];
} Channel;

typedef struct {
	/* ports */
	const LV2_Atom_Sequence* p_control;
	LV2_Atom_Sequence*       p_notify;
	float*                   p_latency;

	/* Config */
	uint32_t n_chn;
	float    rate;
	uint32_t fftlen;
	uint32_t firlen;
	uint32_t parsiz;  // fftlen / 2
	uint32_t firlat;  // firlen / 2
	uint32_t n_segm;  // firlen / fftlen
	uint32_t latency; //parsiz + firlat

	float interp_th; // parsiz / 1e6
	float interp_nm; // 1 / parsiz

	/* Meter Config */
	uint32_t mtr_holdtme;
	uint32_t mtr_fpp;
	float    mtr_falloff;

	/* FFT Buffers and Plan */
	fftwf_complex*  freq_data;
	fftwf_complex*  freq_sum;
	fftwf_complex** freq_fir;
	fftwf_plan      plan_r2c;
	fftwf_plan      plan_c2r;

	/* atom-forge, UI communication */
	LV2_URID_Map*        map;
	ProtLV2URIs          uris;
	LV2_Atom_Forge       forge;
	LV2_Atom_Forge_Frame frame;

	/* GUI state */
	bool  ui_active;
	bool  send_state_to_ui;
	float ui_scale;
	bool  link; // TODO save?

	/* Per channel state */
	Channel channel[MAX_CHANNELS];
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
channel_cleanup (Channel* c)
{
	fftwf_free (c->time_data);
	free (c->buf_dly);
	free (c->buf_src);
	free (c->buf_out);
}

static bool
channel_init (FFTiProc const* self, Channel* c)
{
	c->angle            = 0;
	c->overlap          = 0;
	c->mtr_reset_delay  = 0;
	c->mtr_holdcnt[0]   = 0;
	c->mtr_holdcnt[1]   = 0;
	c->mtr_peak[0]      = 0;
	c->mtr_peak[1]      = 0;
	c->mtr_diff[0]      = 1;
	c->mtr_diff[1]      = 1;
	c->mtr_momentary[0] = 0;
	c->mtr_momentary[1] = 0;

	sin_cos (c->angle, &c->sa, &c->ca);

	c->buf_dly   = (float*)calloc (self->latency, sizeof (float));
	c->buf_src   = (float*)calloc (self->firlen, sizeof (float));
	c->buf_out   = (float*)calloc (self->parsiz, sizeof (float));
	c->time_data = (float*)fftwf_malloc (self->fftlen * sizeof (float));

	return (c->buf_dly && c->buf_src && c->buf_out && c->time_data);
}

static void
channel_reset (FFTiProc const* self, Channel* c)
{
	memset (c->buf_src, 0, self->firlen * sizeof (float));
	memset (c->buf_out, 0, self->parsiz * sizeof (float));
	memset (c->time_data, 0, self->fftlen * sizeof (float));
	c->overlap = 0;
	c->offset  = 0;
}

static void
cleanup (LV2_Handle instance)
{
	FFTiProc* self = (FFTiProc*)instance;

	for (uint32_t chn = 0; chn < MAX_CHANNELS; ++chn) {
		channel_cleanup (&self->channel[chn]);
	}

	fftwf_free (self->freq_data);
	fftwf_free (self->freq_sum);

	if (self->freq_fir) {
		for (uint32_t i = 0; i < self->n_segm; ++i) {
			fftwf_free (self->freq_fir[i]);
		}
	}
	free (self->freq_fir);

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

	if (!strcmp (descriptor->URI, PHASEROTATE_URI)) {
		self->n_chn = 1;
	} else if (!strcmp (descriptor->URI, PHASEROTATE_URI "#stereo")) {
		self->n_chn = 2;
	} else {
		free (self);
		return NULL;
	}

	const LV2_Options_Option* options = NULL;

	for (int i = 0; features[i]; ++i) {
		if (!strcmp (features[i]->URI, LV2_URID__map)) {
			self->map = (LV2_URID_Map*)features[i]->data;
		} else if (!strcmp (features[i]->URI, LV2_OPTIONS__options)) {
			options = (LV2_Options_Option*)features[i]->data;
		}
	}

	if (!self->map) {
		fprintf (stderr, "dpl.lv2 error: Host does not support urid:map\n");
		free (self);
		return NULL;
	}

	lv2_atom_forge_init (&self->forge, self->map);
	map_prot_uris (self->map, &self->uris);

	if (options) {
		LV2_URID atom_Float = self->map->map (self->map->handle, LV2_ATOM__Float);
		LV2_URID ui_scale   = self->map->map (self->map->handle, "http://lv2plug.in/ns/extensions/ui#scaleFactor");
		for (const LV2_Options_Option* o = options; o->key; ++o) {
			if (o->context == LV2_OPTIONS_INSTANCE && o->key == ui_scale && o->type == atom_Float) {
				float ui_scale = *(const float*)o->value;
				if (ui_scale < 1.0) {
					ui_scale = 1.0;
				}
				if (ui_scale > 2.0) {
					ui_scale = 2.0;
				}
				self->ui_scale = ui_scale;
			}
		}
	}

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

	self->rate      = rate;
	self->parsiz    = self->fftlen / 2;
	self->firlat    = self->firlen / 2;
	self->n_segm    = self->firlen / self->parsiz;
	self->interp_th = self->parsiz * 1e-6f;
	self->interp_nm = 1.f / self->parsiz;
	self->latency   = self->parsiz + self->firlat;

	self->ui_active = false;
	self->ui_scale  = 1.0;
	self->link      = false;

	self->mtr_holdtme = 0.5 * rate + 0.5f;
	self->mtr_falloff = 0;
	self->mtr_fpp     = 0;

	self->freq_data = (fftwf_complex*)fftwf_malloc ((self->parsiz + 1) * sizeof (fftwf_complex));
	self->freq_sum  = (fftwf_complex*)fftwf_malloc ((self->parsiz + 1) * sizeof (fftwf_complex));

	bool c_ok = true;
	for (uint32_t chn = 0; chn < self->n_chn; ++chn) {
		c_ok &= channel_init (self, &self->channel[chn]);
	}

	if (!c_ok || !self->freq_data || !self->freq_sum) {
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
	float* time_data = self->channel[0].time_data;
	self->plan_r2c   = fftwf_plan_dft_r2c_1d (self->fftlen, time_data, self->freq_data, FFTW_ESTIMATE);
	self->plan_c2r   = fftwf_plan_dft_c2r_1d (self->fftlen, self->freq_data, time_data, FFTW_ESTIMATE);

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
	memset (&time_data[self->parsiz], 0, sizeof (float) * self->parsiz);

	for (uint32_t s = 0; s < self->n_segm; ++s) {
		for (uint32_t i = 0; i < self->parsiz; ++i) {
			time_data[i] = norm * fir[s * self->parsiz + i];
		}
		fftwf_execute_dft_r2c (self->plan_r2c, time_data, self->freq_fir[s]);
	}

	fftwf_free (freq_hilbert);
	fftwf_free (fir);

	return (LV2_Handle)self;
}

static void
connect_port (LV2_Handle instance,
              uint32_t   port,
              void*      data)
{
	FFTiProc* self = (FFTiProc*)instance;

	switch ((PortIndex)port) {
		case PROT_ATOM_CONTROL:
			self->p_control = (const LV2_Atom_Sequence*)data;
			return;
		case PROT_ATOM_NOTIFY:
			self->p_notify = (LV2_Atom_Sequence*)data;
			return;
		case PROT_LATENCY:
			self->p_latency = (float*)data;
			return;
		default:
			break;
	}

	int chn = (port - PROT_ANGLE0) / 3;
	if (chn < 0 || chn >= MAX_CHANNELS) {
		return;
	}

	Channel* c = &self->channel[chn];
	switch (PROT_ANGLE0 + ((port - PROT_ANGLE0) % 3)) {
		case PROT_INPUT0:
			c->p_in = (float*)data;
			break;
		case PROT_OUTPUT0:
			c->p_out = (float*)data;
			break;
		case PROT_ANGLE0:
			c->p_angle = (float*)data;
			break;
		default:
			break;
	}
}

static float
meter_proc (FFTiProc const* self, Channel* c, float meter_peak, int m)
{
	if (!isfinite (meter_peak)) {
		meter_peak = 0;
	}
	if (meter_peak > c->mtr_peak[m]) {
		c->mtr_peak[m] = meter_peak;
	}
	if (meter_peak > c->mtr_momentary[m]) {
		c->mtr_momentary[m] = meter_peak;
		c->mtr_holdcnt[m]   = self->mtr_holdtme;
	} else if (c->mtr_holdcnt[m] > 0) {
		c->mtr_holdcnt[m] -= self->mtr_fpp;
	} else {
		c->mtr_momentary[m] *= self->mtr_falloff;
		c->mtr_momentary[m] += 1e-20f;
	}
	return meter_peak;
}

static float
meter_run (FFTiProc const* self,
           Channel*        c,
           const float*    d,
           uint32_t        n_samples,
           int             m)
{
	float meter_peak = 0;
	for (uint32_t i = 0; i < n_samples; ++i) {
		float v = fabsf (d[i]);
		if (v > meter_peak) {
			meter_peak = v;
		}
	}
	return meter_proc (self, c, meter_peak, m);
}

static void
meter_reset_peaks (Channel* c)
{
	c->mtr_peak[0] = c->mtr_peak[1] = 0;
	c->mtr_diff[0] = c->mtr_diff[1] = 1;
	c->mtr_momentary[0] = c->mtr_momentary[1] = 0;
}

static void
meter_delayed_reset (Channel* c, uint32_t n_samples, uint32_t reset)
{
	if (c->mtr_reset_delay > 0) {
		c->mtr_diff[0]      = 1;
		c->mtr_diff[1]      = 1;
		c->mtr_momentary[1] = 0;
		c->mtr_reset_delay -= n_samples;
	}
	if (reset > 0) {
		c->mtr_reset_delay = reset + n_samples;
	}
}

static void
activate (LV2_Handle instance)
{
	FFTiProc* self = (FFTiProc*)instance;
	for (uint32_t chn = 0; chn < self->n_chn; ++chn) {
		channel_reset (self, &self->channel[chn]);
		meter_reset_peaks (&self->channel[chn]);
		self->channel[chn].mtr_reset_delay = self->latency;
	}
}

static void
tx_state (FFTiProc* self)
{
	LV2_Atom_Forge_Frame frame;
	lv2_atom_forge_frame_time (&self->forge, 0);
	x_forge_object (&self->forge, &frame, 1, self->uris.state);

	lv2_atom_forge_property_head (&self->forge, self->uris.s_uiscale, 0);
	lv2_atom_forge_float (&self->forge, self->ui_scale);

	lv2_atom_forge_property_head (&self->forge, self->uris.s_link, 0);
	lv2_atom_forge_bool (&self->forge, self->link);

	lv2_atom_forge_pop (&self->forge, &frame);
}

static void
process_channel (FFTiProc* self, uint32_t chn, uint32_t n_samples)
{
	assert (chn < MAX_CHANNELS);
	Channel* c = &self->channel[chn];

	/* localize constants */
	const uint32_t parsiz  = self->parsiz;
	const uint32_t firlen  = self->firlen;
	const uint32_t firlat  = self->firlat;
	const uint32_t n_segm  = self->n_segm;
	const uint32_t latency = self->latency;

	/* localize FFT buffer pointers */
	fftwf_complex*              freq_data = self->freq_data;
	fftwf_complex*              freq_sum  = self->freq_sum;
	fftwf_complex* const* const freq_fir  = self->freq_fir;

	/* per channel settings/data */
	float*   time_data = c->time_data;
	uint32_t offset    = c->offset;
	int32_t  overlap   = c->overlap;
	float    angle     = c->angle;
	uint32_t remain    = n_samples;

	float* iobuf        = c->p_out;
	float  target_angle = *c->p_angle / -360.f;

	if (target_angle < -.5f) {
		target_angle = -.5f;
	}
	if (target_angle > 0.5f) {
		target_angle = 0.5f;
	}

	/* meter delayed input */

	float lvl_in = 0;
	if (n_samples < latency) {
		lvl_in = meter_run (self, c, c->buf_dly, n_samples, 0);
#if 0
		memmove (c->buf_dly, &c->buf_dly[n_samples], n_samples * sizeof (float));
#else
		memcpy (c->buf_dly, &c->buf_dly[n_samples], n_samples * sizeof (float));
#endif
		memcpy (&c->buf_dly[latency - n_samples], iobuf, n_samples * sizeof (float));
	} else {
		float meter_peak = 0;
		for (uint32_t i = 0; i < latency; ++i) {
			float v = fabsf (c->buf_dly[i]);
			if (v > meter_peak) {
				meter_peak = v;
			}
		}
		uint32_t remain = n_samples - latency;
		for (uint32_t i = 0; i < remain; ++i) {
			float v = fabsf (iobuf[i]);
			if (v > meter_peak) {
				meter_peak = v;
			}
		}
		memcpy (c->buf_dly, &iobuf[remain], latency * sizeof (float));
		lvl_in = meter_proc (self, c, meter_peak, 0);
	}

	meter_delayed_reset (c, n_samples, target_angle != angle ? latency : 0);

	/* process output buffer in-place */

	while (remain > 0) {
		uint32_t ns = parsiz - offset;
		if (ns > remain) {
			ns = remain;
		}

		/* copy latent buffers */
		memcpy (&c->buf_src[offset + overlap], iobuf, sizeof (float) * ns);
		memcpy (iobuf, &c->buf_out[offset], sizeof (float) * ns);

		iobuf  += ns;
		offset += ns;
		remain -= ns;

		if (offset == parsiz) {
			offset = 0;

			/* copy end/overlap of prev. iFFT */
			memcpy (c->buf_out, &time_data[parsiz], sizeof (float) * parsiz);

			/* clear convolution buffers */
			memset (&time_data[parsiz], 0, sizeof (float) * parsiz);
			memset (freq_sum, 0, sizeof (fftwf_complex) * (parsiz + 1));

			int32_t olp = overlap;
			for (uint32_t s = 0; s < n_segm; ++s) {
				memcpy (time_data, &c->buf_src[olp], sizeof (float) * parsiz);
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
				c->buf_out[i] += time_data[i];
			}

			/* delayed input (firlen / 2) */
			int32_t off = overlap - firlat;
			if (off < 0) {
				off += firlen;
			}

			float* in  = &c->buf_src[off];
			float* out = c->buf_out;

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
					sin_cos (angle, &c->sa, &c->ca);
				}
			} else {
				const float sa = c->sa;
				const float ca = c->ca;
#pragma GCC ivdep
				for (uint32_t i = 0; i < parsiz; ++i) {
					out[i] = ca * in[i] + sa * out[i];
				}
			}

			overlap = (overlap + parsiz) % firlen;
		}
	}
	/* copy back state */
	c->overlap = overlap;
	c->offset  = offset;
	c->angle   = angle;

	/* Meter Output */
	float lvl_out = meter_run (self, c, c->p_out, n_samples, 1);

	float lvl_diff = 1.0;
	if (c->mtr_momentary[0] > 0.001f && c->mtr_momentary[1] > 0.001f) {
		lvl_diff = c->mtr_momentary[1] / c->mtr_momentary[0];
		if (lvl_diff < c->mtr_diff[0]) {
			c->mtr_diff[0] = lvl_diff;
		}
		if (lvl_diff > c->mtr_diff[1]) {
			c->mtr_diff[1] = lvl_diff;
		}
	}

	if (self->ui_active) {
		LV2_Atom_Forge_Frame frame;
		lv2_atom_forge_frame_time (&self->forge, 0);
		x_forge_object (&self->forge, &frame, 1, self->uris.levels);

		lv2_atom_forge_property_head (&self->forge, self->uris.l_channel, 0);
		lv2_atom_forge_int (&self->forge, chn);

		lv2_atom_forge_property_head (&self->forge, self->uris.l_in_cur, 0);
		lv2_atom_forge_float (&self->forge, lvl_in);
		lv2_atom_forge_property_head (&self->forge, self->uris.l_in_mom, 0);
		lv2_atom_forge_float (&self->forge, c->mtr_momentary[0]);
		lv2_atom_forge_property_head (&self->forge, self->uris.l_in_peak, 0);
		lv2_atom_forge_float (&self->forge, c->mtr_peak[0]);

		lv2_atom_forge_property_head (&self->forge, self->uris.l_out_cur, 0);
		lv2_atom_forge_float (&self->forge, lvl_out);
		lv2_atom_forge_property_head (&self->forge, self->uris.l_out_mom, 0);
		lv2_atom_forge_float (&self->forge, c->mtr_momentary[1]);
		lv2_atom_forge_property_head (&self->forge, self->uris.l_out_peak, 0);
		lv2_atom_forge_float (&self->forge, c->mtr_peak[1]);

		lv2_atom_forge_property_head (&self->forge, self->uris.l_diff_cur, 0);
		lv2_atom_forge_float (&self->forge, lvl_diff);
		lv2_atom_forge_property_head (&self->forge, self->uris.l_diff_min, 0);
		lv2_atom_forge_float (&self->forge, c->mtr_diff[0]);
		lv2_atom_forge_property_head (&self->forge, self->uris.l_diff_max, 0);
		lv2_atom_forge_float (&self->forge, c->mtr_diff[1]);

		lv2_atom_forge_pop (&self->forge, &frame);
	}
}

static void
run (LV2_Handle instance, uint32_t n_samples)
{
	FFTiProc* self = (FFTiProc*)instance;

	/* copy/forward no-inplace buffers */
	for (uint32_t chn = 0; chn < self->n_chn; ++chn) {
		Channel* c = &self->channel[chn];
		if (c->p_in != c->p_out) {
			memcpy (c->p_out, c->p_in, sizeof (float) * n_samples);
		}
	}

	/* announce latency */
	*self->p_latency = self->latency;

	if (!self->p_control || !self->p_notify) {
		/* latency measurement callback */
		return;
	}

	/* prepare forge buffer and initialize atom-sequence */
	const uint32_t capacity = self->p_notify->atom.size;
	lv2_atom_forge_set_buffer (&self->forge, (uint8_t*)self->p_notify, capacity);
	lv2_atom_forge_sequence_head (&self->forge, &self->frame, 0);

	/* process messages from GUI; */
	if (self->p_control) {
		LV2_Atom_Event* ev = lv2_atom_sequence_begin (&(self->p_control)->body);
		while (!lv2_atom_sequence_is_end (&(self->p_control)->body, (self->p_control)->atom.size, ev)) {
			if (ev->body.type == self->uris.atom_Blank || ev->body.type == self->uris.atom_Object) {
				const LV2_Atom_Object* obj = (LV2_Atom_Object*)&ev->body;
				if (obj->body.otype == self->uris.ui_off) {
					self->ui_active = false;
				} else if (obj->body.otype == self->uris.ui_on) {
					self->ui_active        = true;
					self->send_state_to_ui = true;
				} else if (obj->body.otype == self->uris.reset_peaks) {
					for (uint32_t chn = 0; chn < self->n_chn; ++chn) {
						meter_reset_peaks (&self->channel[chn]);
					}
				} else if (obj->body.otype == self->uris.state) {
					const LV2_Atom* v = NULL;
					lv2_atom_object_get (obj, self->uris.s_uiscale, &v, 0);
					if (v) {
						self->ui_scale = ((LV2_Atom_Float*)v)->body;
					}
					v = NULL;
					lv2_atom_object_get (obj, self->uris.s_link, &v, 0);
					if (v) {
						self->link = ((LV2_Atom_Bool*)v)->body;
					}
				}
			}
			ev = lv2_atom_sequence_next (ev);
		}
	}

	/* Calculate meter falloff constants */
	if (self->mtr_fpp != n_samples) {
		const float fall  = 15.0f;
		const float tme   = (float)n_samples / self->rate;     // period time in seconds
		self->mtr_falloff = powf (10.0f, -0.05f * fall * tme); // per period fallback multiplier
		self->mtr_fpp     = n_samples;
	}

	/* Process each Channel */
	for (uint32_t chn = 0; chn < self->n_chn; ++chn) {
		process_channel (self, chn, n_samples);
	}

	if (self->ui_active && self->send_state_to_ui) {
		self->send_state_to_ui = false;
		tx_state (self);
	}

	/* close off atom-sequence */
	lv2_atom_forge_pop (&self->forge, &self->frame);
}

static const void*
extension_data (const char* uri)
{
	return NULL;
}

static LV2_Descriptor descriptor = {
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
			descriptor.URI = PHASEROTATE_URI;
			return &descriptor;
		case 1:
			descriptor.URI = PHASEROTATE_URI "#stereo";
			return &descriptor;
		default:
			return NULL;
	}
}
