phaserotate.lv2
===============

phaserotate.lv2 is an audio phase rotation plugin.

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
