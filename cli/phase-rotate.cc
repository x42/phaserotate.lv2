/*
 * Copyright (C) 2021 Robin Gareus <robin@gareus.org>
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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE // needed for M_PI
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <inttypes.h>
#include <limits>
#include <map>
#include <thread>
#include <vector>

#include <fftw3.h>
#include <sndfile.h>

#include "dsp_peak_calc.h"

#define SUBSAMPLE 2
#define MAXSAMPLE (180 * SUBSAMPLE)

class SinCosLut
{
public:
	SinCosLut ()
	{
		float mp = 2.f * M_PI / SUBSAMPLE / -360.0;
		for (int i = 0; i < MAXSAMPLE; ++i) {
#ifdef __APPLE__
			_sin[i] = sinf (mp * i);
			_cos[i] = cosf (mp * i);
#else
			sincosf (mp * i, &_sin[i], &_cos[i]);
#endif
		}
	}

	void
	sincos (int a, float* s, float* c) const
	{
#if 0
		static float mp = 2.f * M_PI / SUBSAMPLE / -360.0;
		sincosf (a * mp, s, c);
#else
		*s = _sin[a];
		*c = _cos[a];
#endif
	}

private:
	float _sin[MAXSAMPLE];
	float _cos[MAXSAMPLE];
};

SinCosLut scl;

static float
coeff_to_dB (float coeff)
{
	if (coeff < 1e-15) {
		return -std::numeric_limits<float>::infinity ();
	}
	return 20.0f * log10f (coeff);
}

static inline float
calc_peak (float* buf, uint32_t n_samples, float pk)
{
#ifdef _HAVE_DSP_COMPUTE_PEAK
	return dsp_compute_peak (buf, n_samples, pk);
#else
	for (uint32_t i = 0; i < n_samples; ++i) {
		pk = std::max (pk, fabsf (buf[i]));
	}
#endif
	return pk;
}

static inline float
calc_rotated_peak (float const* b0, float const* b1, uint32_t n_samples, float pk, int angle)
{
	float sa, ca;
	scl.sincos (angle, &sa, &ca);

#ifdef _HAVE_DSP_COMPUTE_PEAK
	float x[n_samples];

#pragma GCC ivdep
	for (uint32_t i = 0; i < n_samples; ++i) {
		x[i] = ca * b0[i] + sa * b1[i];
	}
	return calc_peak (x, n_samples, pk);

#else

	for (uint32_t i = 0; i < n_samples; ++i) {
		const float x = ca * b0[i] + sa * b1[i];
		pk = std::max (pk, fabsf (x));
	}
	return pk;
#endif
}

/* ****************************************************************************/

class PhaseRotateProc
{
public:
	PhaseRotateProc (uint32_t blksiz)
		: _fftlen (blksiz * 2)
		, _parsiz (blksiz)
		, _firlen (blksiz / 2)
		, _time_data (0)
		, _freq_data (0)
		, _plan_r2c (0)
		, _plan_c2r (0)
	{
		_time_data = (float*)fftwf_malloc (_fftlen * sizeof (float));
		_freq_data = (fftwf_complex*)fftwf_malloc ((_parsiz + 1) * sizeof (fftwf_complex));
		_ffir_data = (fftwf_complex*)fftwf_malloc ((_parsiz + 1) * sizeof (fftwf_complex));
		_plan_r2c  = fftwf_plan_dft_r2c_1d (_fftlen, _time_data, _freq_data, FFTW_ESTIMATE);
		_plan_c2r  = fftwf_plan_dft_c2r_1d (_fftlen, _freq_data, _time_data, FFTW_ESTIMATE);
		_norm      = 0.5 / _parsiz;

		/* generate hilbert FIR */
		float         fir[_fftlen];
		fftwf_complex freq_fir[_firlen + 1];

		for (uint32_t i = 0; i <= _firlen; ++i) {
			const float re = i & 1 ? -1 : 1;
			freq_fir[i][0] = 0;
			freq_fir[i][1] = re;
		}

		fftwf_plan plan_fir_c2r = fftwf_plan_dft_c2r_1d (_parsiz, freq_fir, fir, FFTW_ESTIMATE);
		fftwf_execute_dft_c2r (plan_fir_c2r, freq_fir, fir);
		fftwf_destroy_plan (plan_fir_c2r);

		const double flen = 1.0 / _parsiz;
		for (uint32_t i = 0; i < _parsiz; ++i) {
			fir[i] *= _norm * (1 - cos (2.0 * M_PI * i * flen));
		}

		memset (fir + _parsiz, 0, _parsiz * sizeof (float));
		fftwf_execute_dft_r2c (_plan_r2c, fir, _ffir_data);
	}

	~PhaseRotateProc ()
	{
		fftwf_free (_time_data);
		fftwf_free (_freq_data);
		fftwf_destroy_plan (_plan_r2c);
		fftwf_destroy_plan (_plan_c2r);
	}

	uint32_t
	parsiz () const
	{
		return _parsiz;
	}

	void
	hilbert (float const* in, float* out, float* o_out)
	{
		const uint32_t parsiz = _parsiz;

		/* copy end/overlap of prev. iFFT */
		memcpy (out, o_out, sizeof (float) * parsiz);

		/* fill fft buffer & zero pad */
		memcpy (_time_data, in + parsiz, sizeof (float) * parsiz);
		memset (_time_data + parsiz, 0, sizeof (float) * parsiz);

		fftwf_execute_dft_r2c (_plan_r2c, _time_data, _freq_data);

		/* FIR convolution, 90deg phase shift */
		const float          norm     = _norm;
		const fftwf_complex* freq_fir = _ffir_data;
#pragma GCC ivdep /* do not check for aliasing, buffers do not overlap */
		for (uint32_t i = 0; i <= parsiz; ++i) {
			const float re   = _freq_data[i][0];
			_freq_data[i][0] = norm * (re * freq_fir[i][0] - _freq_data[i][1] * freq_fir[i][1]);
			_freq_data[i][1] = norm * (re * freq_fir[i][1] + _freq_data[i][1] * freq_fir[i][0]);
		}

		fftwf_execute_dft_c2r (_plan_c2r, _freq_data, _time_data);

#pragma GCC ivdep /* do not check for aliasing, buffers do not overlap */
		for (uint32_t i = 0; i < parsiz; ++i) {
			out[i] += _time_data[i];
		}
		memcpy (o_out, _time_data + parsiz, sizeof (float) * parsiz);
	}

	void
	rotate (float const* in, float* out, int angle)
	{
		float          sa, ca;
		const uint32_t parsiz = _parsiz;
		scl.sincos (angle, &sa, &ca);
		const float* fin = &in[_firlen];
#pragma GCC ivdep
		for (uint32_t i = 0; i < parsiz; ++i) {
			out[i] = ca * fin[i] + sa * out[i];
		}
	}

	void
	process (float const* in, float* out, float* o_out, int angle)
	{
		hilbert (in, out, o_out);
		rotate (in, out, angle);
	}

private:
	PhaseRotateProc (PhaseRotateProc const&) = delete;
	uint32_t       _fftlen;
	uint32_t       _parsiz;
	uint32_t       _firlen;
	float          _norm;
	float*         _time_data;
	fftwf_complex* _freq_data;
	fftwf_complex* _ffir_data;
	fftwf_plan     _plan_r2c;
	fftwf_plan     _plan_c2r;
};

/* ****************************************************************************/

typedef std::vector<std::unique_ptr<PhaseRotateProc>> PRPVec;

class PhaseRotate
{
public:
	PhaseRotate (PRPVec&, uint32_t);
	~PhaseRotate ();

	void reset ();
	void ensure_thread_buffers ();

	void analyze (float const* in, int ang_start = 0, int ang_end = 180, int ang_stride = 1, int chn = -1, bool start = false);
	void apply (float* buf, std::vector<int> const& angles);

	int
	n_channels () const
	{
		return _n_chn;
	}

	int
	blksiz () const
	{
		return _proc.front ()->parsiz ();
	}

	float
	peak (int c, int a) const
	{
		if (c < 0 || (uint32_t)c > _n_chn) {
			return peak_all (a);
		}
		if (a < 0) {
			a += MAXSAMPLE;
		}
		return _peak[c][a % MAXSAMPLE];
	}

	float
	peak_all (int a) const
	{
		float p = 0;
		if (a < 0) {
			a += MAXSAMPLE;
		}
		a = a % MAXSAMPLE;
		for (uint32_t c = 0; c < _n_chn; ++c) {
			p = std::max (p, _peak[c][a]);
		}
		return p;
	}

private:
	void thr_process (float const* in, int ang_start = 0, int ang_end = 180, int ang_stride = 1, int chn = -1, bool start = false);
	void thr_apply (float const* buf, int c, int a);

	PRPVec&  _proc;
	uint32_t _n_chn;
	size_t   _n_threads;
	float**  _buf_out;
	float**  _buf_old;
	float**  _buf_olp;
	float**  _peak;
};

PhaseRotate::PhaseRotate (PRPVec& p, uint32_t n_chn)
	: _proc (p)
	, _n_chn (n_chn)
	, _n_threads (_proc.size ())
{
	uint32_t parsiz = _proc.front ()->parsiz ();

	_buf_old = (float**)calloc (n_chn, sizeof (float*));
	_buf_olp = (float**)calloc (n_chn, sizeof (float*));
	_peak    = (float**)calloc (n_chn, sizeof (float*));

	_buf_out = (float**)calloc (_n_threads, sizeof (float*));
	for (size_t t = 0; t < _n_threads; ++t) {
		_buf_out[t] = (float*)calloc (parsiz, sizeof (float));
	}

	for (uint32_t c = 0; c < n_chn; ++c) {
		_buf_old[c] = (float*)calloc (parsiz, sizeof (float));    // input per channel
		_buf_olp[c] = (float*)calloc (parsiz, sizeof (float));    // overlap per channel
		_peak[c]    = (float*)calloc (MAXSAMPLE, sizeof (float)); // peak per channel per angle
	}
}

PhaseRotate::~PhaseRotate ()
{
	for (size_t t = 0; t < _n_threads; ++t) {
		free (_buf_out[t]);
	}

	for (uint32_t c = 0; c < _n_chn; ++c) {
		free (_buf_old[c]);
		free (_buf_olp[c]);
		free (_peak[c]);
	}

	free (_buf_out);
	free (_buf_old);
	free (_buf_olp);
	free (_peak);
}

void
PhaseRotate::reset ()
{
	uint32_t parsiz = _proc.front ()->parsiz ();
	for (uint32_t c = 0; c < _n_chn; ++c) {
		memset (_buf_old[c], 0, parsiz * sizeof (float));
		memset (_buf_olp[c], 0, parsiz * sizeof (float));
		for (uint32_t a = 0; a < MAXSAMPLE; ++a) {
			_peak[c][a] = 0;
		}
	}
}

void
PhaseRotate::ensure_thread_buffers ()
{
	if (_n_threads >= _proc.size ()) {
		return;
	}
	for (size_t t = 0; t < _n_threads; ++t) {
		free (_buf_out[t]);
	}
	free (_buf_out);

	uint32_t parsiz = _proc.front ()->parsiz ();

	_n_threads = _proc.size ();
	_buf_out   = (float**)calloc (_n_threads, sizeof (float*));
	for (size_t t = 0; t < _n_threads; ++t) {
		_buf_out[t] = (float*)calloc (parsiz, sizeof (float));
	}
}

void
PhaseRotate::thr_process (float const* p_in, int ang_start, int ang_end, int ang_stride, int c, bool start)
{
	uint32_t parsiz = _proc[c]->parsiz ();
	uint32_t fftlen = 2 * parsiz;
	uint32_t firlen = parsiz / 2;

	float tdc[fftlen];
	float hil[parsiz];

	memcpy (tdc, _buf_old[c], parsiz * sizeof (float));
	for (uint32_t i = 0; i < parsiz; ++i) {
		tdc[i + parsiz] = p_in[c + i * _n_chn];
	}

	/* remember overlap */
	memcpy (_buf_old[c], &tdc[parsiz], parsiz * sizeof (float));

	/* apply hilbert FIR */
	_proc[c]->hilbert (tdc, hil, _buf_olp[c]);

	int angle = ang_start;
	while (angle <= ang_end) {
		int a = (angle + MAXSAMPLE) % MAXSAMPLE;
		/* calc peak / angle */
		if (angle == 0) {
			_peak[c][a] = calc_peak (_buf_old[c], parsiz, _peak[c][a]);
		} else {
			memcpy (_buf_out[c], hil, parsiz * sizeof (float));

			if (start) {
				_peak[c][a] = calc_rotated_peak (&tdc[firlen], &_buf_out[c][firlen], firlen, _peak[c][a], a);
			} else {
				_peak[c][a] = calc_rotated_peak (&tdc[firlen], _buf_out[c], parsiz, _peak[c][a], a);
			}
		}
		angle += ang_stride;
		if (angle >= ang_end) {
			break;
		}
	}
}

void
PhaseRotate::analyze (float const* p_in, int ang_start, int ang_end, int ang_stride, int chn, bool start)
{
	uint32_t c0 = (chn < 0) ? 0 : chn;
	uint32_t c1 = (chn < 0) ? _n_chn : chn + 1;

	std::vector<std::thread> threads;
	for (uint32_t c = c0; c < c1; ++c) {
		threads.push_back (std::thread (&PhaseRotate::thr_process, this, p_in, ang_start, ang_end, ang_stride, c, start));
	}
	for (auto& th : threads) {
		th.join ();
	}
}

void
PhaseRotate::thr_apply (float const* buf, int c, int a)
{
	uint32_t parsiz = _proc[c]->parsiz ();
	uint32_t fftlen = 2 * parsiz;

	float tdc[fftlen];

	/* copy overlap */
	memcpy (tdc, _buf_old[c], parsiz * sizeof (float));
	/* de-interleave channel */
	for (uint32_t i = 0; i < parsiz; ++i) {
		tdc[i + parsiz] = buf[c + i * _n_chn];
	}
	/* remember overlap */
	memcpy (_buf_old[c], &tdc[parsiz], parsiz * sizeof (float));

	a = (a + MAXSAMPLE) % MAXSAMPLE;
	_proc[c]->process (tdc, _buf_out[c], _buf_olp[c], a);
}

void
PhaseRotate::apply (float* buf, std::vector<int> const& angles)
{
	std::vector<std::thread> threads;
	for (uint32_t c = 0; c < _n_chn; ++c) {
		threads.push_back (std::thread (&PhaseRotate::thr_apply, this, buf, c, angles[c]));
	}
	for (auto& th : threads) {
		th.join ();
	}

	/* interleave */
	const uint32_t parsiz = _proc.front ()->parsiz ();
	for (uint32_t c = 0; c < _n_chn; ++c) {
		for (uint32_t i = 0; i < parsiz; ++i) {
			buf[c + i * _n_chn] = _buf_out[c][i];
		}
	}
}

/* ****************************************************************************/

static void
usage ()
{
	// help2man compatible format (standard GNU help-text)
	printf ("phase-rotate - Audio File Phase Rotation Util.\n\n");
	printf ("Usage: phase-rotate [ OPTIONS ] <file> [out-file]\n\n");

	/* **** "---------|---------|---------|---------|---------|---------|---------|---------|" */
	printf ("Options:\n"
	        "  -a, --angle <n>[,<n>]*     specify phase angle to apply\n"
	        "  -f, --fftlen <num>         process-block size, freq. resolution\n"
	        "  -h, --help                 display this help and exit\n"
	        "  -l, --link-channels        use downmixed mono peak for analysis\n"
	        "  -s, --stride <num>         analysis step-size\n"
	        "  -v, --verbose              show processing information\n"
	        "  -V, --version              print version information and exit\n"
	        "\n");

	printf ("\n"
	        "This utility analyzes the given audio file to find a phase-rotation\n"
	        "angle that results in minimal digital-peak, while retaining overall\n"
	        "sound and loudness.\n"
	        "\n"
	        "If both input and output file are given, the analysis results applied, and\n"
	        "a new file with optimized phase is written. Otherwise the analysis results\n"
	        "are only printed to standard output.\n"
	        "\n"
	        "Analysis is performed in two steps, first a coarse analysis is performed,\n"
	        "calculating peak for angles distanced `stride' degrees apart. Then local\n"
	        "minimums are explored in a second step.\n"
	        "\n"
	        "Verbose analysis allows to plot the digital peak vs phase-rotation.\n"
	        "The output is in gnuplot(1) data file format.\n"
	        "\n"
	        "If the -a option is specified, no analysis is performed but the given,\n"
	        "phase-angle(s) are directly applied. This requires both input and output\n"
	        "files to be given. If a single angle is given it is applied to all channels\n"
	        "of the file. Otherwise one has to specify the same number of phase-angles as\n"
	        "there are channels in the file.\n"
	        "\n");

	printf ("\n"
	        "Examples:\n"
	        "phase-rotate -l my-music.wav out-file.wav\n\n"
	        "phase-rotate -vv -s 3 my-music.wav\n\n"
	        "phase-rotate -a 10,20 in.wav out.wav\n\n");

	printf ("Report bugs to <https://github.com/x42/phaserotate.lv2/issues>\n"
	        "Website: <https://github.com/x42/phaserotate.lv2/>\n");
	::exit (EXIT_SUCCESS);
}

static void
copy_metadata (SNDFILE* infile, SNDFILE* outfile)
{
	SF_CUES           cues;
	SF_BROADCAST_INFO binfo;

	memset (&cues, 0, sizeof (cues));
	memset (&binfo, 0, sizeof (binfo));

	for (int k = SF_STR_FIRST; k <= SF_STR_LAST; ++k) {
		const char* str = sf_get_string (infile, k);
		if (str != NULL) {
			sf_set_string (outfile, k, str);
		}
	}

	if (sf_command (infile, SFC_GET_CUE, &cues, sizeof (cues)) == SF_TRUE)
		sf_command (outfile, SFC_SET_CUE, &cues, sizeof (cues));

	if (sf_command (infile, SFC_GET_BROADCAST_INFO, &binfo, sizeof (binfo)) == SF_TRUE) {
		sf_command (outfile, SFC_SET_BROADCAST_INFO, &binfo, sizeof (binfo));
	}
}

static void
analyze_file (PhaseRotate& pr, SNDFILE* infile, float* buf, int ang_start, int ang_end, int ang_stride, int chn = -1)
{
	int n_channels = pr.n_channels ();
	int blksiz     = pr.blksiz ();

	bool start = true;
	do {
		sf_count_t n = sf_readf_float (infile, buf, blksiz);
		if (n <= 0) {
			break;
		}
		if (n < blksiz) {
			/* silence pad */
			memset (&buf[n_channels * n], 0, sizeof (float) * n_channels * (blksiz - n));
		}
		pr.analyze (buf, ang_start, ang_end, ang_stride, chn, start);
		start = false;
	} while (1);

	memset (buf, 0, blksiz * n_channels * sizeof (float));
	pr.analyze (buf, ang_start, ang_end, ang_stride, chn, false);
}

int
main (int argc, char** argv)
{
	SF_INFO      nfo;
	SNDFILE*     infile     = NULL;
	SNDFILE*     outfile    = NULL;
	float*       buf        = NULL;
	char*        angles_opt = NULL;
	int          stride     = 12 * SUBSAMPLE;
	int          verbose    = 0;
	FILE*        verbose_fd = stdout;
	bool         find_min   = true;
	bool         link_chn   = false;
	unsigned int blksiz     = 0;

	std::vector<int> angles;

	// TODO pick stride to optimize processor usage

	const char* optstring = "a:f:hls:Vv";

	/* clang-format off */
	const struct option longopts[] = {
		{ "angle",         required_argument, 0, 'a' },
		{ "fftlen",        required_argument, 0, 'f' },
		{ "stride",        required_argument, 0, 's' },
		{ "help",          no_argument,       0, 'h' },
		{ "link-channels", no_argument,       0, 'l' },
		{ "version",       no_argument,       0, 'V' },
		{ "verbose",       no_argument,       0, 'v' },
	};
	/* clang-format on */

	int c = 0;
	while (EOF != (c = getopt_long (argc, argv,
	                                optstring, longopts, (int*)0))) {
		switch (c) {
			case 'a':
				angles_opt = optarg;
				break;

			case 'f':
				blksiz = atoi (optarg);
				break;

			case 'h':
				usage ();
				break;

			case 'l':
				link_chn = true;
				break;

			case 's':
				stride = atoi (optarg);
				break;

			case 'V':
				printf ("phase-rotate version %s\n\n", VERSION);
				printf ("Copyright (C) GPL 2021 Robin Gareus <robin@gareus.org>\n");
				exit (EXIT_SUCCESS);
				break;

			case 'v':
				++verbose;
				break;

			default:
				fprintf (stderr, "Error: unrecognized option. See --help for usage information.\n");
				::exit (EXIT_FAILURE);
				break;
		}
	}

	if (optind + 1 > argc) {
		fprintf (stderr, "Error: Missing parameter. See --help for usage information.\n");
		::exit (EXIT_FAILURE);
	}

	if (stride < 1 || stride > 45 * SUBSAMPLE || (MAXSAMPLE % stride) != 0) {
		fprintf (stderr, "Error: 180 deg is not evenly dividable by given stride.\n");
		::exit (EXIT_FAILURE);
	}

	if (blksiz != 0 && (blksiz < 1024 || blksiz > 32768)) {
		fprintf (stderr, "Error: fft-len is out of bounds; valid range 1024..32768\n");
		::exit (EXIT_FAILURE);
	}

	if (angles_opt && optind + 2 > argc) {
		fprintf (stderr, "Error: -a, --angle option requires an output file to be given.\n");
		::exit (EXIT_FAILURE);
	}

	memset (&nfo, 0, sizeof (SF_INFO));

	if ((infile = sf_open (argv[optind], SFM_READ, &nfo)) == 0) {
		fprintf (stderr, "Cannot open '%s' for reading: ", argv[optind]);
		fputs (sf_strerror (NULL), stderr);
		::exit (EXIT_FAILURE);
	}

	if (find_min && !nfo.seekable) {
		fprintf (stderr, "File '%s' is not seekable. ", argv[optind]);
		::exit (EXIT_FAILURE);
	}

	if (optind + 1 < argc) {
		if ((outfile = sf_open (argv[optind + 1], SFM_WRITE, &nfo)) == 0) {
			fprintf (stderr, "Cannot open '%s' for writing: ", argv[optind + 1]);
			fputs (sf_strerror (NULL), stderr);
			::exit (EXIT_FAILURE);
		}
	}

	if (verbose > 1) {
		verbose_fd = stderr;
	}

	if (verbose > 2) {
		char strbuffer[65536];
		sf_command (infile, SFC_GET_LOG_INFO, strbuffer, 65536);
		fputs (strbuffer, verbose_fd);
	} else if (verbose) {
		fprintf (verbose_fd, "Input File      : %s\n", argv[optind]);
		fprintf (verbose_fd, "Sample Rate     : %d Hz\n", nfo.samplerate);
		fprintf (verbose_fd, "Channels        : %d\n", nfo.channels);
	}

	if (angles_opt) {
		find_min  = false;
		char* ang = angles_opt;
		char* saveptr;
		while ((ang = strtok_r (angles_opt, ",", &saveptr)) != NULL) {
			char*  ep;
			double a = strtod (ang, &ep);
			if (*ep != '\0' || a < -180 || a > 180) {
				fprintf (stderr, "Error: Invalid angle speficied, value needs to be -180 .. +180.\n");
				::exit (EXIT_FAILURE);
			}
			angles_opt = NULL;
			angles.push_back (round (a * (float)SUBSAMPLE));
		}
		if (angles.size () == 1) {
			while (angles.size () < (size_t)nfo.channels) {
				angles.push_back (angles.back ());
			}
		}
		if (angles.size () < (size_t)nfo.channels) {
			fprintf (stderr, "Error: file has more channels than angles were specified.\n");
			::exit (EXIT_FAILURE);
		}
		if (verbose) {
			fprintf (verbose_fd, "# Apply phase-shift\n");
			for (int c = 0; c < nfo.channels; ++c) {
				fprintf (verbose_fd, "Channel: %2d Phase: %5.2f deg\n", c + 1, angles[c] / (float)SUBSAMPLE);
			}
		}
	}

	if (blksiz == 0 || blksiz > 32768) {
		blksiz = nfo.samplerate / 8;
	}

	unsigned power_of_two;
	for (power_of_two = 1; 1U << power_of_two < blksiz; ++power_of_two);
	blksiz = std::min (32768, std::max (1024, 1 << power_of_two));

	if (verbose > 1) {
		fprintf (verbose_fd, "Process block-size %d\n", blksiz);
	}

	buf = (float*)malloc (blksiz * nfo.channels * sizeof (float));

	if (!buf) {
		fprintf (stderr, "Out of memory\n");
		::exit (EXIT_FAILURE);
	}

	{
		PRPVec prp;
		int    n_threads = nfo.channels;

		for (int i = 0; i < n_threads; ++i) {
			auto p = std::unique_ptr<PhaseRotateProc> (new PhaseRotateProc (blksiz));
			prp.push_back (std::move (p));
		}

		PhaseRotate pr (prp, nfo.channels);

		if (find_min) {
			if (verbose > 1) {
				fprintf (verbose_fd, "Analyzing using %d process threads, stride = %d\n", n_threads, stride);
			}

			analyze_file (pr, infile, buf, 0, MAXSAMPLE, stride);

			/*
			 * Consider writing gnuplot file
			 * ```
			 * #!/usr/bin/gnuplot -persist
			 * set terminal qt
			 * set xlabel "Phase Rotation [deg]"
			 * set ylabel "Peak [dBFS]"
			 * set datafile missing "-inf"
			 * plot '-' u 1:3 t "chn-1", '' u 1:4 t "chn-2"
			 *  ... # [data, terminate with 'e\n']
			 * e
			 * ```
			 */

			if (verbose > 1) {
				printf ("# Angle mono-peak");
				for (int c = 0; c < nfo.channels; ++c) {
					printf (" chn-%d", c + 1);
				}
				printf ("\n");
				for (uint32_t a = 0; a < MAXSAMPLE; a += stride) {
					printf ("%.2f %.4f", a / (float)SUBSAMPLE, coeff_to_dB (pr.peak_all (a)));
					for (int c = 0; c < nfo.channels; ++c) {
						printf (" %.4f", coeff_to_dB (pr.peak (c, a)));
					}
					printf ("\n");
				}
			}

			std::map<int, std::vector<int>> mins;
			angles.clear ();

			int   min_angle[(int)nfo.channels];
			float p_min[nfo.channels];
			float r_zro[nfo.channels];
			float r_min[nfo.channels];

			for (int c = 0; c < nfo.channels; ++c) {
				float c_min = std::numeric_limits<float>::infinity ();
				float c_max = 0;
				float range;

				r_zro[c] = pr.peak (c, 0);

				for (uint32_t a = 0; a < MAXSAMPLE; a += stride) {
					c_min = std::min (c_min, pr.peak (link_chn ? -1 : c, a));
					c_max = std::max (c_max, pr.peak (link_chn ? -1 : c, a));
				}

				range = c_max - c_min;
				if (range == 0) {
					mins[0].push_back (c);
					continue;
				}
				if (stride > 1) {
					range *= .07;
					p_min[c] = std::numeric_limits<float>::infinity ();
				} else {
					range    = 0;
					p_min[c] = c_min;
				}

				for (uint32_t a = 0; a < MAXSAMPLE; a += stride) {
					float p = pr.peak (link_chn ? -1 : c, a);
					if (p <= c_min + range) {
						mins[a].push_back (c);
						if (verbose > 1) {
							fprintf (verbose_fd, "Consider min: %f (< %f) chn: %d @ %.2f deg\n", p, c_min + range, c, a / (float)SUBSAMPLE);
						}
					}
				}
			}

			if (stride == 1) {
				for (auto& mp : mins) {
					for (auto& cn : mp.second) {
						min_angle[cn] = mp.first;
						r_min[cn]     = pr.peak (cn, mp.first);
					}
				}
			} else {
				/* search local min */
				int stride_2 = (stride + 1) / 2;

				// TODO: binary search, but use all CPU cores..
				for (auto& mp : mins) {
					if (0 != sf_seek (infile, 0, SEEK_SET)) {
						fprintf (stderr, "Failed to rewind input file\n");
						::exit (EXIT_FAILURE);
					}
					pr.reset ();

					int ma = mp.first;

					analyze_file (pr, infile, buf, ma - stride_2, ma + stride_2 + 1, 1, mp.second.size () > 1 ? -1 : mp.second.front ());

					for (auto& cn : mp.second) {
						for (int a = ma - stride_2; a < ma + stride_2 + 1; ++a) {
							float p = pr.peak (link_chn ? -1 : cn, a);
							if (p <= p_min[cn]) {
								p_min[cn]     = p;
								r_min[cn]     = pr.peak (cn, a);
								min_angle[cn] = (a + MAXSAMPLE) % MAXSAMPLE;
							}

							if (verbose > 1) {
								printf ("%.2f %.4f",
								        ((a + MAXSAMPLE) % MAXSAMPLE) / (float)SUBSAMPLE,
								        coeff_to_dB (pr.peak_all (a)));
								for (int c = 0; c < nfo.channels; ++c) {
									printf (" %.4f", coeff_to_dB (pr.peak (c, a)));
								}
								printf ("\n");
							}
						}
					}
				}
			}

			/* collect results */
			float avg_rotate = 0;
			int   avg_count  = 0;
			for (int c = 0; c < nfo.channels; ++c) {
				if (p_min[c] != std::numeric_limits<float>::infinity ()) {
					avg_rotate += min_angle[c];
					++avg_count;
				}
			}
			avg_rotate /= avg_count;
			float avg_dist = MAXSAMPLE / (float)avg_count;

			/* minimize channel phase distance */
			for (int c = 0; c < nfo.channels; ++c) {
				if (p_min[c] == std::numeric_limits<float>::infinity ()) {
					angles.push_back (0);
				} else {
					if (min_angle[c] > 90 * SUBSAMPLE && fabsf (min_angle[c] - avg_rotate) > avg_dist) {
						min_angle[c] -= MAXSAMPLE;
					} else if (avg_rotate > 90 * SUBSAMPLE) {
						min_angle[c] -= MAXSAMPLE;
					}
					angles.push_back (min_angle[c]);
				}
			}

			/* print results */
			if (!outfile || verbose) {
				fprintf (verbose_fd, "# Result -- Minimize digital peak\n");
				for (int c = 0; c < nfo.channels; ++c) {
					if (p_min[c] == std::numeric_limits<float>::infinity ()) {
						fprintf (verbose_fd, "Channel: %2d Phase:   0 deg # cannot find min.\n", c + 1);
					} else {
						fprintf (verbose_fd, "Channel: %2d Phase: %5.2f deg", c + 1, min_angle[c] / (float)SUBSAMPLE);
						if (min_angle[c] != 0) {
							fprintf (verbose_fd, ", gain: %5.2f dB (att. %4.2f to %4.2f dBFS)",
							         coeff_to_dB (r_zro[c]) - coeff_to_dB (r_min[c]),
							         coeff_to_dB (r_zro[c]), coeff_to_dB (r_min[c]));
						}
						fprintf (verbose_fd, "\n");
					}
				}
			}
		} // find min

		if (outfile) {
			copy_metadata (infile, outfile);
			int n_channels = nfo.channels;

			if (find_min) {
				if (0 != sf_seek (infile, 0, SEEK_SET)) {
					fprintf (stderr, "Failed to rewind input file\n");
					::exit (EXIT_FAILURE);
				}

				pr.reset ();
			}

			const uint32_t latency = blksiz / 2;

			uint32_t pad = 0;
			uint32_t off = latency;
			do {
				sf_count_t n = sf_readf_float (infile, buf, blksiz);
				if (n <= 0) {
					break;
				}

				if (n < latency) {
					pad = blksiz - n;
					/* silence remaining buffer */
					memset (&buf[n_channels * n], 0, sizeof (float) * n_channels * pad);
					pad = latency - n;
					n += pad;
				}

				pr.apply (buf, angles);

				n -= off;

				if (n != sf_writef_float (outfile, &buf[off], n)) {
					fprintf (stderr, "Error writing to output file.\n");
					pad = latency;
					break;
				}
				off = 0;
			} while (1);

			sf_count_t n = (sf_count_t)latency - pad;
			if (n > 0) {
				memset (buf, 0, blksiz * n_channels * sizeof (float));
				pr.apply (buf, angles);

				if (n != sf_writef_float (outfile, buf, n)) {
					fprintf (stderr, "Error writing to output file.\n");
				}
			}
			sf_close (outfile);
		} // write file

	} // scope

	sf_close (infile);
	free (buf);
	fftwf_cleanup ();
	return 0;
}
