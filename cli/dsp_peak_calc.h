/*
 * Copyright (C) 2019-2021 Robin Gareus <robin@gareus.org>
 * Copyright (C) 2020 Ayan Shafqat <ayan.x.shafqat@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdint.h>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#define _HAVE_DSP_COMPUTE_PEAK

static float
dsp_compute_peak (const float* buf, uint32_t n_samples, float current)
{
	float tmp = 0.0f;
	vDSP_maxmgv (buf, (vDSP_Stride)1, &tmp, n_samples);
	return fmaxf (current, tmp);
}

#elif (defined __aarch64__) || (defined __arm__)

#include <arm_acle.h>
#include <arm_neon.h>

#define IS_ALIGNED_TO(ptr, bytes) (((uintptr_t)ptr) % (bytes) == 0)
#define _HAVE_DSP_COMPUTE_PEAK

float
dsp_compute_peak (const float* src, uint32_t nframes, float current)
{
	float32x4_t vc0;

	// Broadcast single value to all elements of the register
	vc0 = vdupq_n_f32 (current);

	// While pointer is not aligned, process one sample at a time
	while (!IS_ALIGNED_TO (src, sizeof (float32x4_t)) && (nframes > 0)) {
		float32x4_t x0;

		x0  = vld1q_dup_f32 (src);
		x0  = vabsq_f32 (x0);
		vc0 = vmaxq_f32 (vc0, x0);

		++src;
		--nframes;
	}

	// SIMD portion with aligned buffers
	do {
		while (nframes >= 8) {
			float32x4_t x0, x1;

			x0 = vld1q_f32 (src + 0);
			x1 = vld1q_f32 (src + 4);

			x0 = vabsq_f32 (x0);
			x1 = vabsq_f32 (x1);

			vc0 = vmaxq_f32 (vc0, x0);
			vc0 = vmaxq_f32 (vc0, x1);

			src += 8;
			nframes -= 8;
		}

		while (nframes >= 4) {
			float32x4_t x0;

			x0 = vld1q_f32 (src);

			x0  = vabsq_f32 (x0);
			vc0 = vmaxq_f32 (vc0, x0);

			src += 4;
			nframes -= 4;
		}

		while (nframes >= 2) {
			float32x2_t x0;
			float32x4_t y0;

			x0 = vld1_f32 (src);        // Load two elements
			x0 = vabs_f32 (x0);         // Compute ABS value
			y0 = vcombine_f32 (x0, x0); // Combine two vectors

			vc0 = vmaxq_f32 (vc0, y0);

			src += 2;
			nframes -= 2;
		}
	} while (0);

	// Do remaining samples one frame at a time
	while (nframes > 0) {
		float32x4_t x0;

		x0  = vld1q_dup_f32 (src);
		x0  = vabsq_f32 (x0);
		vc0 = vmaxq_f32 (vc0, x0);

		++src;
		--nframes;
	}

	// Compute the max in register
	do {
		float32x2_t vlo  = vget_low_f32 (vc0);
		float32x2_t vhi  = vget_high_f32 (vc0);
		float32x2_t max0 = vpmax_f32 (vlo, vhi);
		float32x2_t max1 = vpmax_f32 (max0, max0); // Max is now at max1[0]
		current          = vget_lane_f32 (max1, 0);
	} while (0);

	return current;
}

#undef IS_ALIGNED_TO

#elif (defined __x86_64__) || (defined __i386__) || (defined _M_X64) || (defined _M_IX86) // ARCH_X86

#include <immintrin.h>
#include <xmmintrin.h>

#ifdef __AVX__ // TODO runtime detect AVX
#define _HAVE_DSP_COMPUTE_PEAK

#warning using AVX, this limits available architectures

static float
dsp_compute_peak (const float* src, uint32_t n_samples, float current)
{
	const __m256 ABS_MASK = _mm256_set1_ps (-0.0F);

	// Broadcast the current max value to all elements of the YMM register
	__m256 vmax = _mm256_broadcast_ss (&current);

	// Compute single min/max of unaligned portion until alignment is reached
	while ((((intptr_t)src) % 32 != 0) && n_samples > 0) {
		__m256 vsrc;

		vsrc = _mm256_setzero_ps ();
		vsrc = _mm256_castps128_ps256 (_mm_load_ss (src));
		vsrc = _mm256_andnot_ps (ABS_MASK, vsrc);
		vmax = _mm256_max_ps (vmax, vsrc);

		++src;
		--n_samples;
	}

	// Process the aligned portion 16 samples at a time
	while (n_samples >= 16) {
#ifdef _WIN32
		_mm_prefetch (((char*)src + (16 * sizeof (float))), _mm_hint (0));
#else
		__builtin_prefetch (src + (16 * sizeof (float)), 0, 0);
#endif
		__m256 vsrc1, vsrc2;
		vsrc1 = _mm256_load_ps (src + 0);
		vsrc2 = _mm256_load_ps (src + 8);

		vsrc1 = _mm256_andnot_ps (ABS_MASK, vsrc1);
		vsrc2 = _mm256_andnot_ps (ABS_MASK, vsrc2);

		vmax = _mm256_max_ps (vmax, vsrc1);
		vmax = _mm256_max_ps (vmax, vsrc2);

		src += 16;
		n_samples -= 16;
	}

	// Process the remaining samples 8 at a time
	while (n_samples >= 8) {
		__m256 vsrc;

		vsrc = _mm256_load_ps (src);
		vsrc = _mm256_andnot_ps (ABS_MASK, vsrc);
		vmax = _mm256_max_ps (vmax, vsrc);

		src += 8;
		n_samples -= 8;
	}

	// If there are still some left 4 to 8 samples, process them below
	while (n_samples > 0) {
		__m256 vsrc;

		vsrc = _mm256_setzero_ps ();
		vsrc = _mm256_castps128_ps256 (_mm_load_ss (src));
		vsrc = _mm256_andnot_ps (ABS_MASK, vsrc);
		vmax = _mm256_max_ps (vmax, vsrc);

		++src;
		--n_samples;
	}

	__m256 tmp;
	tmp  = _mm256_shuffle_ps (vmax, vmax, _MM_SHUFFLE (2, 3, 0, 1));
	vmax = _mm256_max_ps (tmp, vmax);
	tmp  = _mm256_shuffle_ps (vmax, vmax, _MM_SHUFFLE (1, 0, 3, 2));
	vmax = _mm256_max_ps (tmp, vmax);
	tmp  = _mm256_permute2f128_ps (vmax, vmax, 1);
	vmax = _mm256_max_ps (tmp, vmax);

	// zero upper 128 bit of 256 bit ymm register to avoid penalties using non-AVX instructions
	_mm256_zeroupper ();

#if defined(__GNUC__) && (__GNUC__ < 5)
	return *((float*)&vmax);
#elif defined(__GNUC__) && (__GNUC__ < 8)
	return vmax[0];
#else
	return _mm256_cvtss_f32 (vmax);
#endif
}

#elif defined __SSE2__
#define _HAVE_DSP_COMPUTE_PEAK

static float
dsp_compute_peak (const float* src, uint32_t n_samples, float current)
{
	const __m128 ABS_MASK = _mm_set1_ps (-0.0F);

	__m128 vmax;
	__m128 temp;

	vmax = _mm_set1_ps (current);

	// Compute single max of unaligned portion until alignment is reached
	while (((intptr_t)src) % 16 != 0 && n_samples > 0) {
		temp = _mm_set1_ps (*src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		++src;
		--n_samples;
	}

	// use 64 byte prefetch for quadruple quads
	while (n_samples >= 16) {
#ifdef _WIN32
		_mm_prefetch (((char*)src + (16 * sizeof (float))), _mm_hint (0));
#else
		__builtin_prefetch (src + (16 * sizeof (float)), 0, 0);
#endif
		temp = _mm_load_ps (src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		src += 4;
		temp = _mm_load_ps (src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		src += 4;
		temp = _mm_load_ps (src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		src += 4;
		temp = _mm_load_ps (src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		src += 4;
		n_samples -= 16;
	}

	// temp through aligned buffers
	while (n_samples >= 4) {
		temp = _mm_load_ps (src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		src += 4;
		n_samples -= 4;
	}

	// temp through the rest < 4 samples
	while (n_samples > 0) {
		temp = _mm_set1_ps (*src);
		temp = _mm_andnot_ps (ABS_MASK, temp);
		vmax = _mm_max_ps (vmax, temp);
		++src;
		--n_samples;
	}

	temp = _mm_shuffle_ps (vmax, vmax, _MM_SHUFFLE (2, 3, 0, 1));
	vmax = _mm_max_ps (temp, vmax);
	temp = _mm_shuffle_ps (vmax, vmax, _MM_SHUFFLE (1, 0, 3, 2));
	vmax = _mm_max_ps (temp, vmax);

#if defined(__GNUC__) && (__GNUC__ < 5)
	return *((float*)&vmax);
#else
	return vmax[0];
#endif
}
#endif // SSE

#endif

#ifndef _HAVE_DSP_COMPUTE_PEAK
#warning dsp_compute_peak is not accelerated on this architecture

static float
dsp_compute_peak (const float* const buf, uint32_t n_samples, float current)
{
	for (uint32_t i = 0; i < n_samples; ++i) {
		const float x = fabsf (buf[i]);
		if (x > current) {
			current = x;
		}
	}
	return current;
}
#endif
