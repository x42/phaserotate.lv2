phaserotate.lv2
===============

phaserotate.lv2 is an audio phase rotation plugin.

Phase shift vs. time shift
--------------------------

Time shifting delays a signal by an **absolute** value in seconds (or samples).
While phase shift uses a **relative** value, a fraction of the wave-length.

When phase shifting a signal, different components of the signal are delayed
differently depending on their frequency.

In the following example two sine-waves are phase shifted by 90deg (1/4 of their
wave-length). Note that the time-delay differs: The sine-wave with the lower
frequency (top, red) is delayed more compared to the higher pitched
sine-wave (bottom, violet).

![](https://github.com/x42/phaserotate.lv2/blob/master/img/phase-rotated-sines-90deg.png "90deg phase-rotated sine-waves")

The interesting aspect is that this does not alter the sound of the signal nor
the loudness. However changing the phase vs. frequency relationship between lower
and upper harmonics changes the waveform and can affect where the peak occurs.

For this reason phase rotation is commonly used by radio stations to reduce
the signal peak and make the signal more symmetrical.

Phase rotation circuits are also used during mastering to increase headroom
in order to add gain and further compress the signal.

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
```

Note to packagers: The Makefile honors `PREFIX` and `DESTDIR` variables as well
as `CFLAGS`, `LDFLAGS` and `OPTIMIZATIONS` (additions to `CFLAGS`), also
see the first 10 lines of the Makefile.
