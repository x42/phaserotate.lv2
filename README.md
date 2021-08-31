phaserotate.lv2
===============

phaserotate.lv2 is an audio phase rotation plugin.

Phase shift vs. time shift
--------------------------

Time shifting delays a signal by an **absolute** value in seconds (or samples).
While phase shift uses a **relative** value, a fraction of the wave-length.

When phase shifting a signal, different components of the signal are delayed
differently depending on their frequency.

The following figure shows an oscilloscope display of two sine-waves with different
frequencies that play concurrently (green and blue) and respective shifted
versions of each signal (red, violet).

The left images shows a phase shift by 90deg (1/4 wave-length).
Note that the time-delay differs: The sine-wave with the lower frequency (top, red)
is delayed further compared to the higher pitched sine-wave (bottom, violet).

In the right image both sine-waves are delayed by the same time. Note how this
results in a different phase shift, as both waves differ in their wavelength.

![](https://github.com/x42/phaserotate.lv2/blob/master/img/time-vs-phase-delay.png "Time vs Phase delay")

Due to the periodic nature of the signal, a phase shift of +180 deg is equivalent to
a shift of -180 deg. The signal is *rotated* around the unit circle on the complex plane.

Phase rotation
--------------

The interesting aspect is that phase rotation does not alter the sound of the signal
nor the loudness. However changing the phase vs. frequency relationship between lower
and upper harmonics changes the waveform and can affect where the digital peak occurs.

For this reason phase rotation is commonly used by radio stations to reduce
the signal peak and make the signal more symmetrical. Phase rotation circuits are
also used during mastering to increase headroom. This allows to increase gain and
further compress the signal.

Caveat
------

By nature of the process phase-rotation affects transient response. Similar to a
linear-phase EQ, there is also pre-ringing. Use this effect with caution because
you are about to become a casualty in the *loudness-war* against quality.
You have been warned.

CLI
---

In addition to the plugin, there is a command-line tool that can analyze an
audio file, and find the phase shift which will result in a minimum peak.

This is for offline analysis and offline processing only, and not installed
by default.

Install
-------

Compiling phaserotate.lv2 requires the LV2 SDK, fftw, gnu-make, and a c-compiler.

```bash
git clone git://github.com/x42/phaserotate.lv2.git
cd phaserotate.lv2
make

ln -s `pwd`/build ~/.lv2/phaserotate.lv2
#sudo make install PREFIX=/usr

cd cli
make
#sudo make install PREFIX=/usr
./phase-rotate --help
```

Note to packagers: The Makefile honors `PREFIX` and `DESTDIR` variables as well
as `CFLAGS`, `LDFLAGS` and `OPTIMIZATIONS` (additions to `CFLAGS`), also
see the first 10 lines of the Makefile.
