/*
 * Copyright (C) 2022 Robin Gareus <robin@gareus.org>
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

#ifndef _PHASEROTATE_H
#define _PHASEROTATE_H

#define PHASEROTATE_URI "http://gareus.org/oss/lv2/phaserotate"

#define PROT_URI PHASEROTATE_URI "#"

#ifdef HAVE_LV2_1_18_6
#include <lv2/atom/atom.h>
#include <lv2/atom/forge.h>
#include <lv2/urid/urid.h>
#else
#include <lv2/lv2plug.in/ns/ext/atom/atom.h>
#include <lv2/lv2plug.in/ns/ext/atom/forge.h>
#include <lv2/lv2plug.in/ns/ext/urid/urid.h>
#endif

#ifdef HAVE_LV2_1_8
#define x_forge_object lv2_atom_forge_object
#else
#define x_forge_object lv2_atom_forge_blank
#endif

typedef struct {
	LV2_URID atom_Blank;
	LV2_URID atom_Object;
	LV2_URID atom_Vector;
	LV2_URID atom_Float;
	LV2_URID atom_Int;
	LV2_URID atom_eventTransfer;
	LV2_URID ui_on;
	LV2_URID ui_off;
	LV2_URID reset_peaks;
	LV2_URID state;
	LV2_URID s_uiscale;
	LV2_URID s_link;
	LV2_URID levels;
	LV2_URID l_channel;
	LV2_URID l_in_cur;
	LV2_URID l_in_mom;
	LV2_URID l_in_peak;
	LV2_URID l_out_cur;
	LV2_URID l_out_mom;
	LV2_URID l_out_peak;
	LV2_URID l_diff_cur;
	LV2_URID l_diff_min;
	LV2_URID l_diff_max;
} ProtLV2URIs;

static inline void
map_prot_uris (LV2_URID_Map* map, ProtLV2URIs* uris)
{
	uris->atom_Blank         = map->map (map->handle, LV2_ATOM__Blank);
	uris->atom_Object        = map->map (map->handle, LV2_ATOM__Object);
	uris->atom_Vector        = map->map (map->handle, LV2_ATOM__Vector);
	uris->atom_Float         = map->map (map->handle, LV2_ATOM__Float);
	uris->atom_Int           = map->map (map->handle, LV2_ATOM__Int);
	uris->atom_eventTransfer = map->map (map->handle, LV2_ATOM__eventTransfer);
	uris->ui_on              = map->map (map->handle, PROT_URI "ui_on");
	uris->ui_off             = map->map (map->handle, PROT_URI "ui_off");
	uris->reset_peaks        = map->map (map->handle, PROT_URI "reset_peaks");
	uris->state              = map->map (map->handle, PROT_URI "state");
	uris->s_uiscale          = map->map (map->handle, PROT_URI "uiscale");
	uris->s_link             = map->map (map->handle, PROT_URI "link");
	uris->levels             = map->map (map->handle, PROT_URI "levels");
	uris->l_channel          = map->map (map->handle, PROT_URI "l_channel");
	uris->l_in_cur           = map->map (map->handle, PROT_URI "l_in_cur");
	uris->l_in_mom           = map->map (map->handle, PROT_URI "l_in_mom");
	uris->l_in_peak          = map->map (map->handle, PROT_URI "l_in_peak");
	uris->l_out_cur          = map->map (map->handle, PROT_URI "l_out_cur");
	uris->l_out_mom          = map->map (map->handle, PROT_URI "l_out_mom");
	uris->l_out_peak         = map->map (map->handle, PROT_URI "l_out_peak");
	uris->l_diff_cur         = map->map (map->handle, PROT_URI "l_diff_cur");
	uris->l_diff_min         = map->map (map->handle, PROT_URI "l_diff_min");
	uris->l_diff_max         = map->map (map->handle, PROT_URI "l_diff_max");
}

/* common definitions UI and DSP */

#define MAX_CHANNELS 2

typedef enum {
	PROT_ATOM_CONTROL = 0,
	PROT_ATOM_NOTIFY,
	PROT_LATENCY,

	PROT_ANGLE0,
	PROT_INPUT0,
	PROT_OUTPUT0,

	PROT_ANGLE1,
	PROT_INPUT1,
	PROT_OUTPUT1
} PortIndex;

#endif
