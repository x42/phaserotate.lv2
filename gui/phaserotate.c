/* robtk phaserotate gui
 *
 * Copyright 2022 Robin Gareus <robin@gareus.org>
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
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lv2/lv2plug.in/ns/ext/atom/atom.h"
#include "lv2/lv2plug.in/ns/ext/options/options.h"

#include "../src/phaserotate.h"

#define RTK_USE_HOST_COLORS
#define RTK_URI PHASEROTATE_URI
#define RTK_GUI "#ui"

#ifndef MAX
#define MAX(A, B) ((A) > (B)) ? (A) : (B)
#endif

#ifndef MIN
#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#endif

typedef struct {
	LV2UI_Write_Function write;
	LV2UI_Controller     controller;
	LV2UI_Touch*         touch;

	LV2_Atom_Forge forge;
	LV2_URID_Map*  map;
	ProtLV2URIs    uris;

	uint32_t n_chn;

	PangoFontDescription* font[2];

	RobWidget* rw;   // top-level container
	RobWidget* ctbl; // control element table

	/* Geometry */
	int  meter_width;
	int  meter_height;
	bool meter_size_changed;
	bool delta_size_changed;
	int  metric_height;

	/* Widgets and Drawing Areas */
	RobWidget* m[6];
	RobWidget* spc_m[3];
	RobTkLbl*  lbl_m[3];

	RobTkDial* spn_ctrl[MAX_CHANNELS];
	RobTkLbl*  lbl_ctrl[MAX_CHANNELS];
	RobTkCBtn* btn_link;

	/* metric area and meter patterns */
	cairo_surface_t* ma[4];
	cairo_pattern_t* lpat;
	cairo_pattern_t* mpat;
	cairo_pattern_t* dpat;
	cairo_pattern_t* spat;

	cairo_surface_t* dial_bg;
	cairo_surface_t* sf_nfo;

	/* meter values */
	float lvl_in[3];
	float lvl_out[3];
	float lvl_diff[3];

	float lvr_in[3];
	float lvr_out[3];
	float lvr_diff[3];

	/* UI */
	bool disable_signals;

	float c_txt[4];
	float c_bgr[4];

	const char* nfo;
} phrotUI;

/* ****************************************************************************/

static float c_dlf[4] = { 0.8, 0.8, 0.8, 1.0 };

/* *****************************************************************************
 * Knob faceplates
 */

#define CTL_WIDTH 165
#define CTL_HEIGHT 126
#define CTL_RADIUS 43
#define CTL_CX 82.5
#define CTL_CY 56.5

static void
prepare_faceplates (phrotUI* ui)
{
	cairo_t* cr;
	float    xlp, ylp;

/* clang-format off */
#define INIT_DIAL_SF(VAR, W, H)                                             \
  VAR = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 2 * (W), 2 * (H)); \
  cr  = cairo_create (VAR);                                                 \
  cairo_scale (cr, 2.0, 2.0);                                               \
  CairoSetSouerceRGBA (c_trs);                                              \
  cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);                           \
  cairo_rectangle (cr, 0, 0, W, H);                                         \
  cairo_fill (cr);                                                          \
  cairo_set_operator (cr, CAIRO_OPERATOR_OVER);

#define DIALDOTS(V, XADD, YADD)                                \
  float ang = (-1. * M_PI) + (2 * M_PI) * (V);                 \
  xlp       = CTL_CX + XADD + sinf (ang) * (CTL_RADIUS + 5.0); \
  ylp       = CTL_CY + YADD - cosf (ang) * (CTL_RADIUS + 5.0); \
  cairo_set_line_cap (cr, CAIRO_LINE_CAP_ROUND);               \
  CairoSetSouerceRGBA (c_dlf);                                 \
  cairo_set_line_width (cr, 3.5);                              \
  cairo_move_to (cr, rint (xlp) - .5, rint (ylp) - .5);        \
  cairo_close_path (cr);                                       \
  cairo_stroke (cr);

#define RESPLABLEL(V)                                             \
  {                                                               \
    DIALDOTS (V, 4.5, 15.5)                                       \
    xlp = rint (CTL_CX + 3.5 + sinf (ang) * (CTL_RADIUS + 12.5));  \
    ylp = rint (CTL_CY + 15.5 - cosf (ang) * (CTL_RADIUS + 12.5)); \
  }
	/* clang-format on */

	INIT_DIAL_SF (ui->dial_bg, CTL_WIDTH + 8, CTL_HEIGHT + 20);
	RESPLABLEL (0.00);
	write_text_full (cr, "+/-180", ui->font[0], xlp, ylp, 0, 2, c_dlf);
	RESPLABLEL (0.125);
	RESPLABLEL (0.25);
	write_text_full (cr, "-90", ui->font[0], xlp - 4, ylp, 0, 2, c_dlf);
	RESPLABLEL (0.375);
	RESPLABLEL (0.5);
	write_text_full (cr, "0", ui->font[0], xlp, ylp, 0, 2, c_dlf);
	RESPLABLEL (.625);
	RESPLABLEL (.75);
	write_text_full (cr, "+90", ui->font[0], xlp + 4, ylp, 0, 2, c_dlf);
	RESPLABLEL (.875);
	cairo_destroy (cr);

#undef DIALDOTS
#undef INIT_DIAL_SF
#undef RESPLABLEL
}

/* *****************************************************************************
 * Numeric value display - knob tooltips
 */

static void
display_annotation (phrotUI* ui, RobTkDial* d, cairo_t* cr, const char* txt)
{
	int tw, th;
	cairo_save (cr);
	PangoLayout* pl = pango_cairo_create_layout (cr);
	pango_layout_set_font_description (pl, ui->font[0]);
	pango_layout_set_text (pl, txt, -1);
	pango_layout_get_pixel_size (pl, &tw, &th);
	cairo_translate (cr, d->w_width / 2, d->w_height / 2);
	cairo_translate (cr, -tw / 2.0, -th / 2.0);
	cairo_set_source_rgba (cr, .0, .0, .0, .7);
	rounded_rectangle (cr, -1, -1, tw + 3, th + 1, 3);
	cairo_fill (cr);
	CairoSetSouerceRGBA (c_wht);
	pango_cairo_show_layout (cr, pl);
	g_object_unref (pl);
	cairo_restore (cr);
	cairo_new_path (cr);
}

static void
dial_annotation (RobTkDial* d, cairo_t* cr, void* data)
{
	phrotUI*    ui = (phrotUI*)(data);
	char        txt[16];
	const float val = d->cur;
	snprintf (txt, 16, "%.1f deg", val);
	display_annotation (ui, d, cr, txt);
}

/* *****************************************************************************
 * Meter drawing
 */

#define MA_HEIGHT (ui->metric_height)

static inline float
deflect_dbfs (int w, float db)
{
	return w * (db + 80.f) / 86.f;
}

static float
deflect_meter (int w, float v)
{
	if (v < .0001f) { // < -80dBFS
		return 0;
	}
	if (v > 2.f) { // > +6.02 dBFS
		return w;
	}
	return deflect_dbfs (w, 20 * log10f (v));
}

static inline float
deflect_db (int w, float db)
{
	return w * (db + 12.f) / 24.f;
}

static float
deflect_delta (int w, float v)
{
	if (v < .252f) { // < -12dB
		return 0;
	}
	if (v > 3.98f) { // > +12 dB
		return w;
	}
	return deflect_db (w, 20 * log10f (v));
}

static void
create_meter_surfaces (phrotUI* ui, int w, int h)
{
	cairo_t* cr;
	char     fn[32];

	snprintf (fn, 32, "Mono %.0fpx", 9. * ui->m[0]->widget_scale);
	PangoFontDescription* font = pango_font_description_from_string (fn);

	if (ui->ma[0]) {
		cairo_surface_destroy (ui->ma[0]);
	}
	if (ui->ma[1]) {
		cairo_surface_destroy (ui->ma[1]);
	}

	/* clang-format off */
#define TICK(DB, TB)                                                    \
  {                                                                     \
    const int x = 1 + deflect_dbfs (w - 2, DB);                         \
    cairo_move_to (cr, .5 + x, TB ? 1.5 : MA_HEIGHT - 1.5);             \
    cairo_close_path (cr);                                              \
    cairo_stroke (cr);                                                  \
    char txt[8];                                                        \
    snprintf (txt, 8, "%+.0f", DB);                                     \
    write_text_full (cr, txt, font, x - 1,                              \
                     TB ? 2 : MA_HEIGHT - 1, 0, TB ? 8 : 5, ui->c_txt); \
  }

#define XTICKS(TB)                                                \
  CairoSetSouerceRGBA (ui->c_txt);                                \
  cairo_set_line_width (cr, 1.5);                                 \
  cairo_set_line_cap (cr, CAIRO_LINE_CAP_ROUND);                  \
  TICK (0.f, TB)                                                  \
  TICK (-6.f, TB)                                                 \
  TICK (-18.f, TB)                                                \
  TICK (-30.f, TB)                                                \
  TICK (-40.f, TB)                                                \
  TICK (-50.f, TB)                                                \
  TICK (-60.f, TB)                                                \
  TICK (-70.f, TB)                                                \
  cairo_set_line_width (cr, 1.0);                                 \
  for (int i = -75; i <= -20; i += 5) {                           \
    const int x = 1 + deflect_dbfs (w - 2, i);                    \
    cairo_move_to (cr, .5 + x, TB ? 1.5 : MA_HEIGHT - 1.5);       \
    cairo_close_path (cr);                                        \
    cairo_stroke (cr);                                            \
  }                                                               \
  for (int i = -15; i <= 3; i += 3) {                             \
    const int x = 1 + deflect_dbfs (w - 2, i);                    \
    cairo_move_to (cr, .5 + x, TB ? 1.5 : MA_HEIGHT - 1.5);       \
    cairo_close_path (cr);                                        \
    cairo_stroke (cr);                                            \
  }                                                               \
  write_text_full (cr, "-\u221E", font, 0, TB ? 1 : MA_HEIGHT - 1, 0, TB ? 9 : 6, ui->c_txt);
	/* clang-format on */

	ui->ma[0] = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, w, MA_HEIGHT);
	cr        = cairo_create (ui->ma[0]);
	cairo_rectangle (cr, 0, 0, w, MA_HEIGHT);
	CairoSetSouerceRGBA (ui->c_bgr);
	cairo_fill (cr);
	XTICKS (false);
	cairo_destroy (cr);

	ui->ma[1] = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, w, MA_HEIGHT);
	cr        = cairo_create (ui->ma[1]);
	cairo_rectangle (cr, 0, 0, w, MA_HEIGHT);
	CairoSetSouerceRGBA (ui->c_bgr);
	cairo_fill (cr);
	XTICKS (true);
	cairo_destroy (cr);

	pango_font_description_free (font);
#undef TICK
#undef XTICKS
}

static void
create_diff_metric_area (phrotUI* ui, int w, int h)
{
	cairo_t* cr;
	char     fn[32];

	snprintf (fn, 32, "Mono %.0fpx", 9. * ui->m[0]->widget_scale);
	PangoFontDescription* font = pango_font_description_from_string (fn);

	if (ui->ma[2]) {
		cairo_surface_destroy (ui->ma[2]);
	}
	if (ui->ma[3]) {
		cairo_surface_destroy (ui->ma[3]);
	}

	/* clang-format off */
#define TICK(DB, TB)                                                    \
  {                                                                     \
    const int x = 1 + deflect_db (w - 2, DB);                           \
    cairo_move_to (cr, .5 + x, TB ? 1.5 : MA_HEIGHT - 1.5);             \
    cairo_close_path (cr);                                              \
    cairo_stroke (cr);                                                  \
    char txt[8];                                                        \
    snprintf (txt, 8, "%+.0f", DB);                                     \
    write_text_full (cr, txt, font, x - 1,                              \
                     TB ? 2 : MA_HEIGHT - 1, 0, TB ? 8 : 5, ui->c_txt); \
  }

#define XTICKS(TB)                                                \
  CairoSetSouerceRGBA (ui->c_txt);                                \
  cairo_set_line_width (cr, 1.5);                                 \
  cairo_set_line_cap (cr, CAIRO_LINE_CAP_ROUND);                  \
  TICK (9.f, TB)                                                 \
  TICK (6.f, TB)                                                  \
  TICK (3.f, TB)                                                  \
  TICK (0.f, TB)                                                  \
  TICK (-3.f, TB)                                                 \
  TICK (-6.f, TB)                                                 \
  TICK (-9.f, TB)                                                \
  cairo_set_line_width (cr, 1.0);                                 \
  for (int i = -11; i <= 11; i += 1) {                            \
    const int x = 1 + deflect_db (w - 2, i);                      \
    cairo_move_to (cr, .5 + x, TB ? 1.5 : MA_HEIGHT - 1.5);       \
    cairo_close_path (cr);                                        \
    cairo_stroke (cr);                                            \
  }
	/* clang-format on */

	ui->ma[2] = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, w, MA_HEIGHT);
	cr        = cairo_create (ui->ma[2]);
	cairo_rectangle (cr, 0, 0, w, MA_HEIGHT);
	CairoSetSouerceRGBA (ui->c_bgr);
	cairo_fill (cr);
	XTICKS (false);
	cairo_destroy (cr);

	ui->ma[3] = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, w, MA_HEIGHT);
	cr        = cairo_create (ui->ma[3]);
	cairo_rectangle (cr, 0, 0, w, MA_HEIGHT);
	CairoSetSouerceRGBA (ui->c_bgr);
	cairo_fill (cr);
	XTICKS (true);
	cairo_destroy (cr);

	pango_font_description_free (font);
#undef TICK
#undef XTICKS
}

static void
create_meter_pattern (phrotUI* ui, int w, int h)
{
	if (ui->lpat) {
		cairo_pattern_destroy (ui->lpat);
	}
	if (ui->mpat) {
		cairo_pattern_destroy (ui->mpat);
	}
	h -= (3 - ui->n_chn) * MA_HEIGHT;

#define PSTOP(DB, R, G, B) \
	cairo_pattern_add_color_stop_rgb (pat, deflect_dbfs (w - 2, DB) / (w - 2.0), R, G, B);

	cairo_pattern_t* pat = cairo_pattern_create_linear (0.0, 0.0, w - 2, 0);
	cairo_pattern_add_color_stop_rgb (pat, .0, .0, .0, .0);
	/* clang-format off */
	PSTOP (-70.0, 0.0, 0.2, 0.5);
	PSTOP (-65.0, 0.0, 0.5, 0.2);
	PSTOP (-18.3, 0.0, 0.7, 0.0);
	PSTOP (-18.0, 0.0, 1.0, 0.0);
	PSTOP ( -9.3, 0.0, 1.0, 0.0);
	PSTOP ( -9.0, 0.7, 0.7, 0.0);
	PSTOP ( -3.3, 0.7, 0.7, 0.0);
	PSTOP ( -3.0, 0.8, 0.5, 0.0);
	PSTOP ( -0.3, 1.0, 0.5, 0.0);
	PSTOP (  0.0, 1.0, 0.0, 0.0);
	/* clang-format on */
	cairo_pattern_add_color_stop_rgb (pat, 1.0, 1.0, .0, .0);

	{
		cairo_pattern_t* shade_pattern = cairo_pattern_create_linear (0.0, 0, 0, h);
		cairo_pattern_add_color_stop_rgba (shade_pattern, 0.0, 0.0, 0.0, 0.0, 0.90);
		cairo_pattern_add_color_stop_rgba (shade_pattern, .26, 0.0, 0.0, 0.0, 0.55);
		cairo_pattern_add_color_stop_rgba (shade_pattern, .40, 1.0, 1.0, 1.0, 0.12);
		cairo_pattern_add_color_stop_rgba (shade_pattern, .53, 0.0, 0.0, 0.0, 0.05);
		cairo_pattern_add_color_stop_rgba (shade_pattern, .74, 0.0, 0.0, 0.0, 0.55);
		cairo_pattern_add_color_stop_rgba (shade_pattern, 1.0, 0.0, 0.0, 0.0, 0.90);

		cairo_surface_t* surface;
		surface     = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, w, h);
		cairo_t* tc = cairo_create (surface);
		cairo_set_source (tc, pat);
		cairo_rectangle (tc, 0, 0, w, h);
		cairo_fill (tc);

		cairo_set_source (tc, shade_pattern);
		cairo_rectangle (tc, 0, 0, w, h);
		cairo_fill (tc);
		cairo_pattern_destroy (shade_pattern);

		ui->lpat = cairo_pattern_create_for_surface (surface);
		cairo_destroy (tc);
		cairo_surface_destroy (surface);
	}
	ui->mpat = pat;

#undef PSTOP
}

static void
create_delta_pattern (phrotUI* ui, int w, int h)
{
	if (ui->dpat) {
		cairo_pattern_destroy (ui->dpat);
	}
	if (ui->spat) {
		cairo_pattern_destroy (ui->spat);
	}
	h -= (3 - ui->n_chn) * MA_HEIGHT;

#define PSTOP(DB, R, G, B) \
	cairo_pattern_add_color_stop_rgb (pat, deflect_db (w - 2, DB) / (w - 2.0), R, G, B);

	cairo_pattern_t* pat = cairo_pattern_create_linear (0.0, 0.0, w - 2, 0);
	cairo_pattern_add_color_stop_rgb (pat, .0, 1.0, .0, .0);
	/* clang-format off */
	PSTOP (-12.0, 1.0, 0.0, 0.0);
	PSTOP ( -6.1, 0.8, 0.5, 0.0);
	PSTOP ( -5.9, 0.8, 0.8, 0.0);
	PSTOP ( -3.0, 0.8, 0.8, 0.0);
	PSTOP ( -2.9, 0.0, 1.0, 0.0);
	PSTOP ( -1.1, 0.0, 1.0, 0.0);
	PSTOP ( -0.9, 0.0, 0.7, 0.0);
	PSTOP (  0.0, 0.0, 0.5, 0.2);
	PSTOP (  0.9, 0.0, 0.7, 0.0);
	PSTOP (  1.1, 0.0, 1.0, 0.0);
	PSTOP (  2.9, 0.0, 1.0, 0.0);
	PSTOP (  3.0, 0.8, 0.8, 0.0);
	PSTOP (  5.9, 0.8, 0.8, 0.0);
	PSTOP (  6.1, 0.8, 0.5, 0.0);
	PSTOP ( 12.0, 1.0, 0.0, 0.0);
	/* clang-format on */
	cairo_pattern_add_color_stop_rgb (pat, 1.0, 1.0, .0, .0);

	{
		cairo_pattern_t* shade_pattern = cairo_pattern_create_linear (0.0, 0, 0, h);
		cairo_pattern_add_color_stop_rgba (shade_pattern, 0.0, 0.0, 0.0, 0.0, 0.2);
		cairo_pattern_add_color_stop_rgba (shade_pattern, 1.0, 1.0, 1.0, 1.0, 0.4);

		cairo_surface_t* surface;
		surface     = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, w, h);
		cairo_t* tc = cairo_create (surface);
		cairo_set_source (tc, pat);
		cairo_rectangle (tc, 0, 0, w, h);
		cairo_fill (tc);

		cairo_set_source (tc, shade_pattern);
		cairo_rectangle (tc, 0, 0, w, h);
		cairo_fill (tc);
		cairo_pattern_destroy (shade_pattern);

		ui->dpat = cairo_pattern_create_for_surface (surface);
		cairo_destroy (tc);
		cairo_surface_destroy (surface);
	}

	cairo_pattern_destroy (pat);

	cairo_pattern_t* shade_pattern = cairo_pattern_create_linear (0.0, 0.0, 0.0, h);
	cairo_pattern_add_color_stop_rgba (shade_pattern, 0.0, 0.0, 0.0, 0.0, 0.90);
	cairo_pattern_add_color_stop_rgba (shade_pattern, .26, 0.0, 0.0, 0.0, 0.55);
	cairo_pattern_add_color_stop_rgba (shade_pattern, .40, 1.0, 1.0, 1.0, 0.12);
	cairo_pattern_add_color_stop_rgba (shade_pattern, .53, 0.0, 0.0, 0.0, 0.05);
	cairo_pattern_add_color_stop_rgba (shade_pattern, .74, 0.0, 0.0, 0.0, 0.55);
	cairo_pattern_add_color_stop_rgba (shade_pattern, 1.0, 0.0, 0.0, 0.0, 0.90);
	ui->spat = shade_pattern;
#undef PSTOP
}

static void
level_expose_event (phrotUI* ui, cairo_t* cr, cairo_rectangle_t* ev, float* lvl, int ma)
{
	int w = ui->meter_width;
	int h = ui->meter_height;
	int y = 0;

	if (ui->meter_size_changed) {
		ui->meter_size_changed = false;
		create_meter_surfaces (ui, w, h);
		create_meter_pattern (ui, w, h);
	}

	cairo_rectangle (cr, ev->x, ev->y, ev->width, ev->height);
	cairo_clip_preserve (cr);
	CairoSetSouerceRGBA (ui->c_bgr);
	cairo_fill (cr);

	/* metric area */
	if (ma & 1) {
		cairo_set_source_surface (cr, ui->ma[0], 0, 0);
		cairo_paint (cr);
	}

	if (ma & 2) {
		cairo_set_source_surface (cr, ui->ma[1], 0, h - MA_HEIGHT);
		cairo_paint (cr);
		h -= MA_HEIGHT;
	}

	if (ma & 1) {
		h -= MA_HEIGHT;
		y += MA_HEIGHT;
	}

	/* meter border/bg */
	cairo_save (cr);

	rounded_rectangle (cr, 0.5, y + 0.5, w - 1, h - 1, 6);
	CairoSetSouerceRGBA (c_gry);
	cairo_set_line_width (cr, 1.0);
	cairo_stroke (cr);

	rounded_rectangle (cr, 1, y + 1, w - 2, h - 2, 6);
	CairoSetSouerceRGBA (c_blk);
	cairo_fill_preserve (cr);
	cairo_clip (cr);

	/* meter levels */

	float px_lvl, px_mom, px_pik;

	px_lvl = deflect_meter (w - 2, lvl[0]);
	px_mom = deflect_meter (w - 2, lvl[1]);
	px_pik = deflect_meter (w - 2, lvl[2]);

	cairo_translate (cr, 1, y);

	if (px_lvl > 2) {
		rounded_rectangle (cr, 0, 1, px_lvl - 1, h - 2, 6);
		cairo_set_source (cr, ui->lpat);
		cairo_fill (cr);
	}

	if (px_mom > 0) {
		rounded_rectangle (cr, -.5 + px_mom, 1, .5, h - 2, 6);
		cairo_set_source (cr, ui->mpat);
		cairo_fill (cr);
	}

	if (px_pik > 0) {
		cairo_set_line_cap (cr, CAIRO_LINE_CAP_BUTT);
		cairo_set_line_width (cr, 2.0);
		cairo_move_to (cr, px_pik, 1);
		cairo_rel_line_to (cr, 0, h - 2);
		cairo_set_source (cr, ui->lpat);
		cairo_stroke_preserve (cr);
		cairo_set_source_rgba (cr, 1, .7, .5, 0.3);
		cairo_stroke (cr);
	}
	cairo_restore (cr);
}

static bool
gain_expose_event (phrotUI* ui, cairo_t* cr, cairo_rectangle_t* ev, float* lvl, int ma)
{
	int w = ui->meter_width;
	int h = ui->meter_height;
	int y = 0;

	if (ui->delta_size_changed) {
		ui->delta_size_changed = false;
		create_diff_metric_area (ui, w, h);
		create_delta_pattern (ui, w, h);
	}

	cairo_rectangle (cr, ev->x, ev->y, ev->width, ev->height);
	cairo_clip_preserve (cr);
	CairoSetSouerceRGBA (ui->c_bgr);
	cairo_fill (cr);

	/* metric area */
	if (ma & 1) {
		cairo_set_source_surface (cr, ui->ma[2], 0, 0);
		cairo_paint (cr);
	}

	if (ma & 2) {
		cairo_set_source_surface (cr, ui->ma[3], 0, h - MA_HEIGHT);
		cairo_paint (cr);
		h -= MA_HEIGHT;
	}

	if (ma & 1) {
		h -= MA_HEIGHT;
		y += MA_HEIGHT;
	}

	int w2 = (w - 2) / 2;

	/* meter border/bg */
	cairo_save (cr);

	rounded_rectangle (cr, 0.5, y + 0.5, w - 1, h - 1, 6);
	CairoSetSouerceRGBA (c_gry);
	cairo_set_line_width (cr, 1.0);
	cairo_stroke (cr);

	rounded_rectangle (cr, 1, y + 1, w - 2, h - 2, 6);
	CairoSetSouerceRGBA (c_blk);
	cairo_fill_preserve (cr);
	cairo_clip (cr);

	/* delta levels */

	const float px_lvl = deflect_delta (w - 2, lvl[0]);
	const float px_min = deflect_delta (w - 2, lvl[1]);
	const float px_max = deflect_delta (w - 2, lvl[2]);

	cairo_translate (cr, 1, y);
	cairo_set_source (cr, ui->dpat);

	if (px_min < w2 - 4 || px_max > w2 + 4) {
		cairo_rectangle (cr, px_min, 1, px_max - px_min, h - 2);
		cairo_fill_preserve (cr);
		cairo_set_source_rgba (cr, 0, 0, 0, .25);
		cairo_fill (cr);
	}

	if (px_lvl > w2 + 2) {
		int r = MIN (6, px_lvl - w2);
		int y = h / 4;
		cairo_new_sub_path (cr);
		cairo_arc (cr, px_lvl - r, y + r, r, M_PI / -2.0, 0);
		cairo_arc (cr, px_lvl - r, h - r - y, r, 0, M_PI / 2.0);
		cairo_line_to (cr, w2, h - y);
		cairo_line_to (cr, w2, y);
		cairo_close_path (cr);
		cairo_fill_preserve (cr);
		cairo_set_source (cr, ui->spat);
		cairo_fill (cr);
	} else if (px_lvl < w2 - 2) {
		int r = MIN (6, w2 - px_lvl);
		int y = h / 4;
		cairo_new_sub_path (cr);
		cairo_line_to (cr, w2, y);
		cairo_line_to (cr, w2, h - y);
		cairo_arc (cr, r + px_lvl, h - r - y, r, M_PI / 2.0, M_PI);
		cairo_arc (cr, r + px_lvl, y + r, r, M_PI, 1.5 * M_PI);
		cairo_close_path (cr);
		cairo_fill_preserve (cr);
		cairo_set_source (cr, ui->spat);
		cairo_fill (cr);
	} else {
		cairo_set_line_cap (cr, CAIRO_LINE_CAP_BUTT);
		cairo_set_line_width (cr, 2.0);
		cairo_move_to (cr, w2, 3);
		cairo_rel_line_to (cr, 0, h - 6);
		cairo_stroke (cr);
	}

	double dash = 1;
	cairo_set_dash (cr, &dash, 1, 1);
	cairo_set_line_cap (cr, CAIRO_LINE_CAP_BUTT);
	cairo_set_line_width (cr, 1.0);
	cairo_move_to (cr, .5 + w2, 1);
	cairo_rel_line_to (cr, 0, h - 2);
	cairo_set_source_rgba (cr, 1, 1, 1, .7);
	cairo_stroke (cr);

	cairo_restore (cr);

	return TRUE;
}

static bool
meter_expose_event (RobWidget* handle, cairo_t* cr, cairo_rectangle_t* ev)
{
	phrotUI* ui = (phrotUI*)GET_HANDLE (handle);
	if (handle == ui->m[0]) {
		level_expose_event (ui, cr, ev, ui->lvl_in, ui->n_chn == 1 ? 3 : 1);
	} else if (handle == ui->m[1]) {
		level_expose_event (ui, cr, ev, ui->lvl_out, ui->n_chn == 1 ? 3 : 1);
	} else if (handle == ui->m[2]) {
		gain_expose_event (ui, cr, ev, ui->lvl_diff, ui->n_chn == 1 ? 3 : 1);
	} else if (handle == ui->m[3]) {
		level_expose_event (ui, cr, ev, ui->lvr_in, 2);
	} else if (handle == ui->m[4]) {
		level_expose_event (ui, cr, ev, ui->lvr_out, 2);
	} else if (handle == ui->m[5]) {
		gain_expose_event (ui, cr, ev, ui->lvr_diff, 2);
	} else {
		assert (0);
	}
	return TRUE;
}

static void
meter_size_request (RobWidget* handle, int* w, int* h)
{
	phrotUI* ui       = (phrotUI*)GET_HANDLE (handle);
	ui->metric_height = ui->m[0]->widget_scale * 12;
	*w                = ui->m[0]->widget_scale * 240;
	*h                = ui->m[0]->widget_scale * (ui->n_chn > 1 ? 11 : 17);
	*h               += MA_HEIGHT * (3 - ui->n_chn);
}

static void
meter_size_allocate (RobWidget* m, int w, int h)
{
	phrotUI* ui = (phrotUI*)GET_HANDLE (m);
	if (m == ui->m[0]) {
		ui->meter_width        = w;
		ui->meter_height       = h;
		ui->meter_size_changed = true;
		ui->delta_size_changed = true;
		if (ui->sf_nfo) {
			cairo_surface_destroy (ui->sf_nfo);
			ui->sf_nfo = NULL;
		}
	}
	assert (w == ui->meter_width);
	assert (h == ui->meter_height);
	robwidget_set_size (m, w, h);
	queue_draw (m);
}

static bool
spacer_expose_event (RobWidget* handle, cairo_t* cr, cairo_rectangle_t* ev)
{
	return TRUE;
}

static void
spacer_size_request (RobWidget* handle, int* w, int* h)
{
	phrotUI* ui = (phrotUI*)GET_HANDLE (handle);
	*w          = ui->m[0]->widget_scale * 2;
	*h          = ui->m[0]->widget_scale * 4;
}

static bool
xtable_expose_event (RobWidget* rw, cairo_t* cr, cairo_rectangle_t* ev)
{
	phrotUI* ui = (phrotUI*)rw->top;
	rcontainer_expose_event (rw, cr, ev);

	if (!ui->nfo) {
		return TRUE;
	}
	if (!ui->sf_nfo) {
		int   tw, th;
		float c[4];
		char  fn[32];
		get_color_from_theme (1, c);
		snprintf (fn, 32, "Sans %.0fpx", 10. * ui->m[0]->widget_scale);
		PangoFontDescription* fd = pango_font_description_from_string (fn);
		get_text_geometry (ui->nfo, fd, &tw, &th);
		/* create_text_surface2 w/backround */
		ui->sf_nfo  = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, th, tw);
		cairo_t* cr = cairo_create (ui->sf_nfo);
		CairoSetSouerceRGBA (c);
		cairo_set_operator (cr, CAIRO_OPERATOR_SOURCE);
		cairo_paint (cr);
		cairo_set_operator (cr, CAIRO_OPERATOR_OVER);
		write_text_full (cr, ui->nfo, fd, 0, tw, 1.5 * M_PI, 9, c_g40);
		cairo_surface_flush (ui->sf_nfo);
		cairo_destroy (cr);
		pango_font_description_free (fd);
	}

	if (ui->sf_nfo) {
		cairo_set_source_surface (cr, ui->sf_nfo, 1, 1);
		cairo_paint (cr);
	}

	return TRUE;
}

/* *****************************************************************************
 * knob callbacks
 */

static bool
cb_spn_ctrl (RobWidget* w, void* handle)
{
	phrotUI* ui = (phrotUI*)handle;

	if (ui->disable_signals) {
		return TRUE;
	}

	if (robtk_cbtn_get_active (ui->btn_link)) {
		if (w == ui->spn_ctrl[0]->rw) {
			const float val = robtk_dial_get_value (ui->spn_ctrl[0]);
			robtk_dial_set_value (ui->spn_ctrl[1], val);
		}
	}

	for (uint32_t i = 0; i < ui->n_chn; ++i) {
		if (w == ui->spn_ctrl[i]->rw) {
			const float val = robtk_dial_get_value (ui->spn_ctrl[i]);
			ui->write (ui->controller, PROT_ANGLE0 + 3 * i, sizeof (float), 0, (const void*)&val);
			break;
		}
	}
	return TRUE;
}

static bool
btn_link (RobWidget* w, void* handle)
{
	phrotUI* ui = (phrotUI*)handle;
	if (robtk_cbtn_get_active (ui->btn_link)) {
		robtk_dial_set_sensitive (ui->spn_ctrl[1], false);
		cb_spn_ctrl (ui->spn_ctrl[0]->rw, handle);
	} else {
		robtk_dial_set_sensitive (ui->spn_ctrl[1], true);
	}
	return TRUE;
}

static RobWidget*
meter_mousedown (RobWidget* handle, RobTkBtnEvent* event)
{
	phrotUI* ui = (phrotUI*)GET_HANDLE (handle);

	uint8_t obj_buf[64];
	lv2_atom_forge_set_buffer (&ui->forge, obj_buf, 64);
	LV2_Atom_Forge_Frame frame;
	lv2_atom_forge_frame_time (&ui->forge, 0);
	LV2_Atom* msg = (LV2_Atom*)x_forge_object (&ui->forge, &frame, 1, ui->uris.reset_peaks);
	lv2_atom_forge_pop (&ui->forge, &frame);
	ui->write (ui->controller, PROT_ATOM_CONTROL, lv2_atom_total_size (msg), ui->uris.atom_eventTransfer, msg);

	return NULL;
}

/* ****************************************************************************/

static RobWidget*
toplevel (phrotUI* ui, void* const top)
{
	/* main widget: layout */
	ui->rw = rob_vbox_new (FALSE, 2);
	robwidget_make_toplevel (ui->rw, top);
	robwidget_toplevel_enable_scaling (ui->rw);

	ui->font[0] = pango_font_description_from_string ("Mono 9px");

	prepare_faceplates (ui);

	/* control knob table */
	ui->ctbl      = rob_table_new (/*rows*/ 1, /*cols*/ 1, FALSE);
	ui->ctbl->top = (void*)ui;

	robwidget_set_expose_event (ui->ctbl, xtable_expose_event);

#define GSP_W(PTR) robtk_dial_widget (PTR)
#define GLB_W(PTR) robtk_lbl_widget (PTR)
#define GBT_W(PTR) robtk_cbtn_widget (PTR)

	for (uint32_t i = 0; i < ui->n_chn; ++i) {
		if (ui->n_chn == 2) {
			ui->lbl_ctrl[i] = robtk_lbl_new (i == 0 ? "Angle Left" : "Angle Right");
		} else {
			ui->lbl_ctrl[i] = robtk_lbl_new ("Angle");
		}

		ui->spn_ctrl[i] = robtk_dial_new_with_size (
		    -180, 180, 0.5,
		    CTL_WIDTH + 8, CTL_HEIGHT + 20, CTL_CX + 4, CTL_CY + 15, CTL_RADIUS);

		robtk_dial_set_value (ui->spn_ctrl[i], 0);
		robtk_dial_set_callback (ui->spn_ctrl[i], cb_spn_ctrl, ui);
		robtk_dial_set_default (ui->spn_ctrl[i], 0);
		robtk_dial_set_scroll_mult (ui->spn_ctrl[i], 10);

		if (ui->touch) {
			robtk_dial_set_touch (ui->spn_ctrl[i], ui->touch->touch, ui->touch->handle, PROT_ANGLE0);
		}

		robtk_dial_set_scaled_surface_scale (ui->spn_ctrl[i], ui->dial_bg, 2.0);
		/* snap at 0 degree */
		robtk_dial_set_detent_default (ui->spn_ctrl[i], true);

		/* use 'dot' for time knobs */
		ui->spn_ctrl[i]->displaymode       = 3; // use dot
		ui->spn_ctrl[i]->threesixty        = true;
		ui->spn_ctrl[i]->with_scroll_accel = false;

		/* show numeric */
		robtk_dial_annotation_callback (ui->spn_ctrl[i], dial_annotation, ui);
	}

	/* link button */
	ui->btn_link = robtk_cbtn_new ("Link", GBT_LED_OFF, false);
	robtk_cbtn_set_alignment (ui->btn_link, 0.5, 0.5);
	robtk_cbtn_set_callback (ui->btn_link, btn_link, ui);
	robtk_cbtn_set_color_checked (ui->btn_link, .45, .45, .45);

	/* meters */
	for (uint32_t i = 0; i < 3 * ui->n_chn; ++i) {
		ui->m[i] = robwidget_new (ui);
		robwidget_set_expose_event (ui->m[i], meter_expose_event);
		robwidget_set_size_request (ui->m[i], meter_size_request);
		robwidget_set_size_allocate (ui->m[i], meter_size_allocate);
		robwidget_set_mousedown (ui->m[i], meter_mousedown);
	}

	ui->meter_size_changed = true;
	ui->delta_size_changed = true;

	ui->lbl_m[0] = robtk_lbl_new ("Input");
	ui->lbl_m[1] = robtk_lbl_new ("Output");
	ui->lbl_m[2] = robtk_lbl_new ("Gain");

	for (uint32_t i = 0; i < 3; ++i) {
		robtk_lbl_set_alignment (ui->lbl_m[i], 1.0, .47);
	}

	for (uint32_t i = 0; i < 3; ++i) {
		ui->spc_m[i] = robwidget_new (ui);
		robwidget_set_size_request (ui->spc_m[i], spacer_size_request);
		robwidget_set_expose_event (ui->spc_m[i], spacer_expose_event);
	}

	if (ui->n_chn == 1) {
		rob_table_attach (ui->ctbl, GSP_W (ui->spn_ctrl[0]), 0, 1, 0, 3, 0, 0, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_ctrl[0]), 0, 1, 3, 4, 0, 0, RTK_SHRINK, RTK_SHRINK);

		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_m[0]), 1, 2, 0, 1, 0, 8, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_m[2]), 1, 2, 1, 2, 0, 8, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_m[1]), 1, 2, 2, 4, 0, 8, RTK_SHRINK, RTK_SHRINK);

		rob_table_attach (ui->ctbl, ui->m[0], 2, 3, 0, 1, 2, 8, RTK_EXANDF, RTK_SHRINK);
		rob_table_attach (ui->ctbl, ui->m[2], 2, 3, 1, 2, 2, 8, RTK_EXANDF, RTK_SHRINK);
		rob_table_attach (ui->ctbl, ui->m[1], 2, 3, 2, 4, 2, 8, RTK_EXANDF, RTK_SHRINK);
	} else {
		/* clang-format off */
		rob_table_attach (ui->ctbl, GSP_W (ui->spn_ctrl[0]), 0, 1,  0, 10, 0, 0, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_ctrl[0]), 0, 1, 10, 11, 0, 0, RTK_SHRINK, RTK_SHRINK);

		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_m[0]), 1, 2, 1, 3, 0, 8, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_m[2]), 1, 2, 4, 6, 0, 8, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_m[1]), 1, 2, 7, 9, 0, 8, RTK_SHRINK, RTK_SHRINK);

		rob_table_attach (ui->ctbl, ui->spc_m[0], 2, 3,  0,  1, 2, 1, RTK_EXANDF, RTK_EXANDF);

		rob_table_attach (ui->ctbl, ui->m[0],     2, 3,  1,  2, 2, 1, RTK_EXANDF, RTK_SHRINK);
		rob_table_attach (ui->ctbl, ui->m[3],     2, 3,  2,  3, 2, 1, RTK_EXANDF, RTK_SHRINK);

		rob_table_attach (ui->ctbl, ui->spc_m[1], 2, 3,  3,  4, 2, 1, RTK_EXANDF, RTK_SHRINK);

		rob_table_attach (ui->ctbl, ui->m[2],     2, 3,  4,  5, 2, 1, RTK_EXANDF, RTK_SHRINK);
		rob_table_attach (ui->ctbl, ui->m[5],     2, 3,  5,  6, 2, 1, RTK_EXANDF, RTK_SHRINK);

		rob_table_attach (ui->ctbl, ui->spc_m[2], 2, 3,  6,  7, 2, 1, RTK_EXANDF, RTK_SHRINK);

		rob_table_attach (ui->ctbl, ui->m[1],     2, 3,  7,  8, 2, 1, RTK_EXANDF, RTK_SHRINK);
		rob_table_attach (ui->ctbl, ui->m[4],     2, 3,  8,  9, 2, 1, RTK_EXANDF, RTK_SHRINK);

		rob_table_attach (ui->ctbl, GSP_W (ui->spn_ctrl[1]), 3, 4,  0, 10, 0, 0, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GLB_W (ui->lbl_ctrl[1]), 3, 4, 10, 11, 0, 0, RTK_SHRINK, RTK_SHRINK);
		rob_table_attach (ui->ctbl, GBT_W (ui->btn_link), 1, 3, 10, 11, 0, 0, RTK_EXANDF, RTK_SHRINK);
		/* clang-format on */
	}

	/* top-level packing */
	rob_vbox_child_pack (ui->rw, ui->ctbl, TRUE, TRUE);
	return ui->rw;
}

static void
gui_cleanup (phrotUI* ui)
{
	for (uint32_t i = 0; i < ui->n_chn; ++i) {
		robtk_dial_destroy (ui->spn_ctrl[i]);
		robtk_lbl_destroy (ui->lbl_ctrl[i]);
	}
	for (uint32_t i = 0; i < 3 * ui->n_chn; ++i) {
		robwidget_destroy (ui->m[i]);
	}
	for (uint32_t i = 0; i < 3; ++i) {
		robtk_lbl_destroy (ui->lbl_m[i]);
	}
	for (uint32_t i = 0; i < 3; ++i) {
		robwidget_destroy (ui->spc_m[i]);
	}

	robtk_cbtn_destroy (ui->btn_link);

	cairo_surface_destroy (ui->dial_bg);

	if (ui->sf_nfo) {
		cairo_surface_destroy (ui->sf_nfo);
	}

	cairo_pattern_destroy (ui->lpat);
	cairo_pattern_destroy (ui->mpat);
	cairo_pattern_destroy (ui->dpat);
	cairo_pattern_destroy (ui->spat);
	cairo_surface_destroy (ui->ma[0]);
	cairo_surface_destroy (ui->ma[1]);
	cairo_surface_destroy (ui->ma[2]);
	cairo_surface_destroy (ui->ma[3]);

	pango_font_description_free (ui->font[0]);

	rob_table_destroy (ui->ctbl);
	rob_box_destroy (ui->rw);
}

/* *****************************************************************************
 * RobTk + LV2
 */

#define LVGL_RESIZEABLE

static enum LVGLResize
plugin_scale_mode (LV2UI_Handle handle)
{
	return LVGL_LAYOUT_TO_FIT;
}

static void
tx_state (phrotUI* ui)
{
	uint8_t obj_buf[1024];
	lv2_atom_forge_set_buffer (&ui->forge, obj_buf, 1024);
	LV2_Atom_Forge_Frame frame;
	lv2_atom_forge_frame_time (&ui->forge, 0);
	LV2_Atom* msg = (LV2_Atom*)x_forge_object (&ui->forge, &frame, 1, ui->uris.state);

	lv2_atom_forge_property_head (&ui->forge, ui->uris.s_uiscale, 0);
	lv2_atom_forge_float (&ui->forge, ui->rw->widget_scale);

	lv2_atom_forge_property_head (&ui->forge, ui->uris.s_link, 0);
	lv2_atom_forge_bool (&ui->forge, robtk_cbtn_get_active (ui->btn_link));

	lv2_atom_forge_pop (&ui->forge, &frame);
	ui->write (ui->controller, PROT_ATOM_CONTROL, lv2_atom_total_size (msg), ui->uris.atom_eventTransfer, msg);
}

static void
ui_enable (LV2UI_Handle handle)
{
	phrotUI* ui = (phrotUI*)handle;

	uint8_t obj_buf[64];
	lv2_atom_forge_set_buffer (&ui->forge, obj_buf, 64);
	LV2_Atom_Forge_Frame frame;
	lv2_atom_forge_frame_time (&ui->forge, 0);
	LV2_Atom* msg = (LV2_Atom*)x_forge_object (&ui->forge, &frame, 1, ui->uris.ui_on);
	lv2_atom_forge_pop (&ui->forge, &frame);
	ui->write (ui->controller, PROT_ATOM_CONTROL, lv2_atom_total_size (msg), ui->uris.atom_eventTransfer, msg);
}

static void
ui_disable (LV2UI_Handle handle)
{
	phrotUI* ui = (phrotUI*)handle;

	tx_state (ui); // too late?

	uint8_t obj_buf[64];
	lv2_atom_forge_set_buffer (&ui->forge, obj_buf, 64);
	LV2_Atom_Forge_Frame frame;
	lv2_atom_forge_frame_time (&ui->forge, 0);
	LV2_Atom* msg = (LV2_Atom*)x_forge_object (&ui->forge, &frame, 1, ui->uris.ui_off);
	lv2_atom_forge_pop (&ui->forge, &frame);
	ui->write (ui->controller, PROT_ATOM_CONTROL, lv2_atom_total_size (msg), ui->uris.atom_eventTransfer, msg);
}

static LV2UI_Handle
instantiate (
    void* const               ui_toplevel,
    const LV2UI_Descriptor*   descriptor,
    const char*               plugin_uri,
    const char*               bundle_path,
    LV2UI_Write_Function      write_function,
    LV2UI_Controller          controller,
    RobWidget**               widget,
    const LV2_Feature* const* features)
{
	phrotUI* ui = (phrotUI*)calloc (1, sizeof (phrotUI));
	if (!ui) {
		return NULL;
	}

	if (!strcmp (plugin_uri, RTK_URI)) {
		ui->n_chn = 1;
	} else if (!strcmp (plugin_uri, RTK_URI "#stereo")) {
		ui->n_chn = 2;
	} else {
		free (ui);
		return NULL;
	}

	const LV2_Options_Option* options = NULL;

	for (int i = 0; features[i]; ++i) {
		if (!strcmp (features[i]->URI, LV2_UI__touch)) {
			ui->touch = (LV2UI_Touch*)features[i]->data;
		} else if (!strcmp (features[i]->URI, LV2_URID__map)) {
			ui->map = (LV2_URID_Map*)features[i]->data;
		} else if (!strcmp (features[i]->URI, LV2_OPTIONS__options)) {
			options = (LV2_Options_Option*)features[i]->data;
		}
	}

	if (!ui->map) {
		fprintf (stderr, "phaserotate.lv2 UI: Host does not support urid:map\n");
		free (ui);
		return NULL;
	}

	map_prot_uris (ui->map, &ui->uris);
	lv2_atom_forge_init (&ui->forge, ui->map);
	interpolate_fg_bg (c_dlf, .2);

	get_color_from_theme (0, ui->c_txt);
	get_color_from_theme (1, ui->c_bgr);

	ui->nfo             = robtk_info (ui_toplevel);
	ui->write           = write_function;
	ui->controller      = controller;
	ui->disable_signals = true;
	*widget             = toplevel (ui, ui_toplevel);
	ui->disable_signals = false;

	if (options) {
		LV2_URID atom_Float = ui->map->map (ui->map->handle, LV2_ATOM__Float);
		LV2_URID ui_scale   = ui->map->map (ui->map->handle, "http://lv2plug.in/ns/extensions/ui#scaleFactor");
		for (const LV2_Options_Option* o = options; o->key; ++o) {
			if (o->context == LV2_OPTIONS_INSTANCE && o->key == ui_scale && o->type == atom_Float) {
				float ui_scale = *(const float*)o->value;
				if (ui_scale < 1.0) {
					ui_scale = 1.0;
				}
				if (ui_scale > 2.0) {
					ui_scale = 2.0;
				}
				robtk_queue_scale_change (ui->rw, ui_scale);
			}
		}
	}
	return ui;
}

static void
cleanup (LV2UI_Handle handle)
{
	phrotUI* ui = (phrotUI*)handle;
	gui_cleanup (ui);
	free (ui);
}

static void
set_meter_level (float* lvl, float c, float m, float p, RobWidget* rw)
{
	bool diff = false;
	diff |= fabsf (lvl[0] - c) > 1e-4;
	diff |= fabsf (lvl[1] - m) > 1e-4;
	diff |= fabsf (lvl[2] - p) > 1e-4;
	lvl[0] = c;
	lvl[1] = m;
	lvl[2] = p;
	if (diff && rw) {
		queue_draw (rw);
	}
}

/* receive information from DSP */
static void
port_event (LV2UI_Handle handle,
            uint32_t     port_index,
            uint32_t     buffer_size,
            uint32_t     format,
            const void*  buffer)
{
	phrotUI* ui = (phrotUI*)handle;

	if (format == ui->uris.atom_eventTransfer && port_index == PROT_ATOM_NOTIFY) {
		LV2_Atom* atom = (LV2_Atom*)buffer;
		if (atom->type != ui->uris.atom_Blank && atom->type != ui->uris.atom_Object) {
			return;
		}
		LV2_Atom_Object* obj = (LV2_Atom_Object*)atom;

		if (obj->body.otype == ui->uris.state) {
			const LV2_Atom* a0 = NULL;
			const LV2_Atom* a1 = NULL;
			if (2 == lv2_atom_object_get (obj,
			                              ui->uris.s_uiscale, &a0,
			                              ui->uris.s_link, &a1,
			                              NULL)) {
				const float sc = ((LV2_Atom_Float*)a0)->body;
				if (sc != ui->rw->widget_scale && sc >= 1.0 && sc <= 2.0) {
					robtk_queue_scale_change (ui->rw, sc);
				}
				const bool link = ((LV2_Atom_Bool*)a1)->body;
				robtk_cbtn_set_active (ui->btn_link, link);
			}
		} else if (obj->body.otype == ui->uris.levels) {
			const LV2_Atom *ac, *a0, *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8;
			ac = a0 = a1 = a2 = a3 = a4 = a5 = a6 = a7 = a8 = NULL;
			if (10 == lv2_atom_object_get (obj,
			                               ui->uris.l_channel, &ac,
			                               ui->uris.l_in_cur, &a0,
			                               ui->uris.l_in_mom, &a1,
			                               ui->uris.l_in_peak, &a2,
			                               ui->uris.l_out_cur, &a3,
			                               ui->uris.l_out_mom, &a4,
			                               ui->uris.l_out_peak, &a5,
			                               ui->uris.l_diff_cur, &a6,
			                               ui->uris.l_diff_min, &a7,
			                               ui->uris.l_diff_max, &a8,
			                               NULL)) {
				int channel = ((LV2_Atom_Int*)ac)->body;
#define AF(val) ((LV2_Atom_Float*)val)->body
				if (channel == 0) {
					set_meter_level (ui->lvl_in, AF(a0), AF(a1), AF(a2), ui->m[0]);
					set_meter_level (ui->lvl_out, AF(a3), AF(a4), AF(a5), ui->m[1]);
					set_meter_level (ui->lvl_diff, AF(a6), AF(a7), AF(a8), ui->m[2]);
				} else {
					set_meter_level (ui->lvr_in, AF(a0), AF(a1), AF(a2), ui->m[3]);
					set_meter_level (ui->lvr_out, AF(a3), AF(a4), AF(a5), ui->m[4]);
					set_meter_level (ui->lvr_diff, AF(a6), AF(a7), AF(a8), ui->m[5]);
				}
#undef AF
			}
		}
		return;
	}

	if (format != 0) {
		return;
	}

	int chn = (port_index - PROT_ANGLE0) / 3;
	if (chn >= 0 && chn < (int)ui->n_chn && 0 == (port_index - PROT_ANGLE0) % 3) {
		const float v       = *(float*)buffer;
		ui->disable_signals = true;
		robtk_dial_set_value (ui->spn_ctrl[chn], v);
		ui->disable_signals = false;
	}
}

static const void*
extension_data (const char* uri)
{
	return NULL;
}
