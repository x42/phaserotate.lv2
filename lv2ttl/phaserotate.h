#define MULTIPLUGIN 1
#define X42_MULTIPLUGIN_NAME "Phaserotate"
#define X42_MULTIPLUGIN_URI "http://gareus.org/oss/lv2/phaserotate"

#include "lv2ttl/phaserotate_mono.h"
#include "lv2ttl/phaserotate_stereo.h"

static const RtkLv2Description _plugins[] = {
	_plugin_mono,
	_plugin_stereo,
};
