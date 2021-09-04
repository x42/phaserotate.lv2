PREFIX ?= /usr/local
bindir = $(PREFIX)/bin
mandir = $(PREFIX)/share/man/man1

PKG_CONFIG ?= pkg-config
OPTIMIZATIONS ?= -ffast-math -fomit-frame-pointer -O3 -fno-finite-math-only -DNDEBUG # -fopt-info-vec-missed -fopt-info -fopt-info-loop-vec-all
CXXFLAGS ?= $(OPTIMIZATIONS) -Wall -g -std=c++11

VERSION?=$(shell (git describe --tags HEAD 2>/dev/null || echo "v0.4") | sed 's/^v//')

###############################################################################

ifeq ($(shell $(PKG_CONFIG) --exists sndfile || echo no), no)
  $(error "http://www.mega-nerd.com/libsndfile/ is required - install libsndfile1-dev")
endif

ifeq ($(shell $(PKG_CONFIG) --exists fftw3f || echo no), no)
  $(error "fftw3f library was not found")
endif

CPPFLAGS+=-DVERSION=\"$(VERSION)\"
CXXFLAGS+=`$(PKG_CONFIG) --cflags sndfile fftw3f` -pthread
LOADLIBES=`$(PKG_CONFIG) --libs sndfile fftw3f` -lm

UNAME=$(shell uname)
MACHINE=$(shell uname -m)
ifeq ($(UNAME),Darwin)
  CXXFLAGS+= -framework Accelerate
else
  ifneq (,$(findstring arm,$(MACHINE)))
    CXXFLAGS+= -mfpu=neon-vfpv4
  endif
  ifeq ($(MACHINE),x86_64)
    CXXFLAGS+= -msse -msse2 -mavx -mfpmath=sse
  endif
  ifeq ($(MACHINE),i386)
    CXXFLAGS+= -msse -msse2 -mfpmath=sse
  endif
endif

all: phase-rotate

man: phase-rotate.1

phase-rotate: phase-rotate.cc

phase-rotate.1: phase-rotate
	help2man -N -n 'Audio File Phase Rotator' -o phase-rotate.1 ./phase-rotate

clean:
	rm -f phase-rotate

install: install-bin install-man

uninstall: uninstall-bin uninstall-man

install-bin: phase-rotate
	install -d $(DESTDIR)$(bindir)
	install -m755 phase-rotate $(DESTDIR)$(bindir)

uninstall-bin:
	rm -f $(DESTDIR)$(bindir)/phase-rotate
	-rmdir $(DESTDIR)$(bindir)

install-man:
	install -d $(DESTDIR)$(mandir)
	install -m644 phase-rotate.1 $(DESTDIR)$(mandir)

uninstall-man:
	rm -f $(DESTDIR)$(mandir)/phase-rotate.1
	-rmdir $(DESTDIR)$(mandir)

.PHONY: all clean install uninstall man install-man install-bin uninstall-man uninstall-bin
