"""Shared chrome for the Eswatini presentation decks.

Theme resolution, the BMGF template wiring, and the layout primitives that
`build_deck.py` (the full technical deck) and `build_update.py` (the short
internal update) both draw with. Importing this module parses `--theme` and
binds the palette, so `from deck_kit import *` gives a caller the right colours
at import time rather than after an init() call.

Themes: `bmgf` (default) builds onto the real Gates Foundation template;
`neutral` builds a standalone deck with the dataviz reference palette.
"""

import argparse
import copy
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

def _up(start, test):
    """Nearest ancestor of `start` (inclusive) satisfying `test`."""
    for d in [start, *start.parents]:
        if test(d):
            return d
    raise FileNotFoundError(f"no ancestor of {start} matched")


PRES = Path(__file__).resolve().parent
# The model repo is whichever ancestor holds experiments/ -- resolved rather
# than counted, so moving presentations/ does not break the scripts again.
REPO = _up(PRES, lambda d: (d / "experiments").is_dir())
EXP = REPO / "experiments"
ROOT = REPO

# The Gates Foundation corporate template. Not a file we author -- it is
# lifted from `Parameter_comparison_HIVsim_EMOD.pptx`, which carries the
# `BMGF Layouts` theme and its 70 slide layouts. Building onto it means the
# master, colour scheme, type scale and logo are inherited rather than
# reconstructed from memory, which is the only way to get brand fidelity right.
BMGF_TEMPLATE = _up(
    PRES,
    lambda d: (d / "Parameter_comparison_HIVsim_EMOD.pptx").is_file()
) / "Parameter_comparison_HIVsim_EMOD.pptx"

# --- Themes -------------------------------------------------------------------
# Colours must stay in step with make_figures.py, since the deck's own chrome
# sits next to charts drawn in the same palette.
#
# The BMGF values are read out of the template's theme1.xml colour scheme
# ("Custom 4"), not recalled:
#   lt1 Parchment #F5F3ED | accent6 Slate #303A44 | dk2 Saffron #E7C700
#   accent1 Orange #F85C02 | accent2 Red #D93027 | accent4 Blue #248AF9
#   accent5 Turquoise #3AC9B1
# Its type scale is Calibri Light for headings, Calibri for body.
THEMES = {
    "neutral": dict(
        blue=0x2A78D6, orange=0xEB6834, aqua=0x1BAF7A, red=0xE34948,
        ink=0x22221F, ink2=0x6B6A63, rule=0xD9D8D1, band=0xF2F1EC,
        dark=0x22221F, dark_sub=0xC3C2B7, dark_kicker=0x1BAF7A,
        head_font="Segoe UI", body_font="Segoe UI",
        template=None, figdir="", suffix="",
    ),
    "bmgf": dict(
        blue=0x248AF9, orange=0xF85C02, aqua=0x3AC9B1, red=0xD93027,
        ink=0x303A44, ink2=0x6E7681, rule=0xC7C4BA, band=0xF5F3ED,
        dark=0x303A44, dark_sub=0xD5D2C8, dark_kicker=0xE7C700,
        head_font="Calibri Light", body_font="Calibri",
        template=BMGF_TEMPLATE, figdir="bmgf", suffix="_bmgf",
    ),
}

_p = argparse.ArgumentParser()
_p.add_argument("--theme", choices=sorted(THEMES), default="bmgf")
THEME = _p.parse_known_args()[0].theme
T = THEMES[THEME]

FIG = PRES / "figures" / T["figdir"]
SUFFIX = T["suffix"]          # callers compose their own output filename

_c = lambda v: RGBColor((v >> 16) & 0xFF, (v >> 8) & 0xFF, v & 0xFF)  # noqa: E731
BLUE, ORANGE, AQUA, RED = (_c(T[k]) for k in ("blue", "orange", "aqua", "red"))
INK, INK2 = _c(T["ink"]), _c(T["ink2"])
RULE, BAND = _c(T["rule"]), _c(T["band"])
DARK, DARK_SUB, DARK_KICKER = (_c(T[k]) for k in
                               ("dark", "dark_sub", "dark_kicker"))
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
# Row tint for 'we come back to this one'. Deliberately pale: the text
# on top stays near-black, so contrast is unaffected.
HILITE = _c(0xE4F4EF) if THEME == "bmgf" else _c(0xE8F3EE)
HEAD_FONT, BODY_FONT = T["head_font"], T["body_font"]

W, H = Inches(13.333), Inches(7.5)
MARGIN = Inches(0.62)
BODY_W = W - 2 * MARGIN


def _layout_by_name(pres, name):
    for lay in pres.slide_masters[0].slide_layouts:
        if lay.name == name:
            return lay
    raise KeyError(f"layout {name!r} not in template")


def _strip_slides(pres):
    """Drop the template's own slides, keeping its masters and layouts."""
    ids = pres.slides._sldIdLst
    for sld_id in list(ids):
        pres.part.drop_rel(sld_id.rId)
        ids.remove(sld_id)


if T["template"] is None:
    prs = Presentation()
    prs.slide_width, prs.slide_height = W, H
    BLANK = prs.slide_layouts[6]
    LAY_TITLE = LAY_DIVIDER = LAY_END = None
else:
    prs = Presentation(str(T["template"]))
    _strip_slides(prs)
    prs.slide_width, prs.slide_height = W, H
    BLANK = _layout_by_name(prs, "Blank slide - white")
    LAY_TITLE = _layout_by_name(prs, "Title Slide - Logo Frame, Parchment")
    LAY_DIVIDER = _layout_by_name(prs, "Section Divider, Slate")
    LAY_END = _layout_by_name(prs, "End Slide - Slate")


# --------------------------------------------------------------------------
# primitives
# --------------------------------------------------------------------------
def textbox(slide, left, top, width, height, runs, align=PP_ALIGN.LEFT,
            anchor=MSO_ANCHOR.TOP, spacing=1.0, font=None):
    """runs = list of (text, size, bold, color) or (text, size, bold, color, space_before)."""
    tb = slide.shapes.add_textbox(left, top, width, height)
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, spec in enumerate(runs):
        text, size, bold, color = spec[:4]
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        if len(spec) > 4 and spec[4]:
            p.space_before = Pt(spec[4])
        r = p.add_run()
        r.text = text
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = color
        r.font.name = font or BODY_FONT
    return tb


def rect(slide, left, top, width, height, color):
    from pptx.enum.shapes import MSO_SHAPE
    s = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, left, top, width, height)
    s.fill.solid()
    s.fill.fore_color.rgb = color
    s.line.fill.background()
    s.shadow.inherit = False
    return s


def slide_header(slide, eyebrow, title, accent=BLUE):
    """Eyebrow (experiment reference) + title + rule. Returns y below the rule."""
    top = Inches(0.42)
    if eyebrow:
        textbox(slide, MARGIN, top, BODY_W, Inches(0.24),
                [(eyebrow.upper(), 11, True, accent)])
        top = top + Inches(0.30)
    tb = textbox(slide, MARGIN, top, BODY_W, Inches(1.0),
                 [(title, 25, True, INK)], spacing=0.95,
                 font=HEAD_FONT)
    # ~70 characters is where 25pt Segoe UI Semibold wraps at this width.
    lines = 1 + len(title) // 70
    top = top + Inches(0.42 * lines)
    rect(slide, MARGIN, top + Inches(0.10), BODY_W, Emu(9525), RULE)
    return top + Inches(0.30)


def takeaway(slide, text, accent=BLUE):
    """The one-line reading of the slide, in a band at the bottom."""
    h = Inches(0.72)
    top = H - MARGIN - h
    rect(slide, MARGIN, top, BODY_W, h, BAND)
    rect(slide, MARGIN, top, Inches(0.055), h, accent)
    textbox(slide, MARGIN + Inches(0.26), top, BODY_W - Inches(0.5), h,
            [(text, 14, False, INK)], anchor=MSO_ANCHOR.MIDDLE, spacing=1.15)
    return top


def picture(slide, path, top, bottom_limit, left=None, width=None):
    """Place an image, scaled to fit the box between `top` and `bottom_limit`."""
    from PIL import Image
    iw, ih = Image.open(path).size
    box_w = width if width is not None else BODY_W
    box_h = bottom_limit - top - Inches(0.10)
    scale = min(box_w / iw, box_h / ih)
    w, h = int(iw * scale), int(ih * scale)
    x = left if left is not None else int(MARGIN + (BODY_W - w) / 2)
    if left is not None and width is not None:
        x = int(left + (width - w) / 2)
    slide.shapes.add_picture(str(path), x, int(top + (box_h - h) / 2), w, h)


def table(slide, rows, top, col_w, left=None, size=11.5,
          head_color=INK, row_h=Inches(0.34), highlight=None,
          tint=None):
    """rows[0] is the header. col_w in Inches units summing to <= BODY_W."""
    n_r, n_c = len(rows), len(rows[0])
    total_w = sum(col_w)
    shape = slide.shapes.add_table(n_r, n_c, MARGIN if left is None else left,
                                   top, total_w, row_h * n_r)
    tbl = shape.table
    for j, cw in enumerate(col_w):
        tbl.columns[j].width = cw
    for i, row in enumerate(rows):
        tbl.rows[i].height = row_h if i else Inches(0.32)
        for j, val in enumerate(row):
            cell = tbl.cell(i, j)
            cell.margin_left = Inches(0.08)
            cell.margin_right = Inches(0.06)
            cell.margin_top = cell.margin_bottom = Inches(0.03)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            cell.fill.solid()
            if i == 0:
                cell.fill.fore_color.rgb = BAND
            elif highlight and (i - 1) in highlight:
                # Rows the deck goes on to cover in detail.
                cell.fill.fore_color.rgb = tint or HILITE
            else:
                cell.fill.fore_color.rgb = WHITE
            tf = cell.text_frame
            tf.word_wrap = True
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT
            r = p.add_run()
            r.text = str(val)
            r.font.size = Pt(size)
            r.font.bold = bool(i == 0 or (highlight and (i - 1) in highlight and j == 0))
            r.font.color.rgb = head_color if i == 0 else INK
            r.font.name = BODY_FONT
    return shape


def footnote(slide, text, top):
    textbox(slide, MARGIN, top, BODY_W, Inches(0.4),
            [(text, 9.5, False, INK2)], spacing=1.25)


def new(eyebrow, title, accent=BLUE):
    s = prs.slides.add_slide(BLANK)
    y = slide_header(s, eyebrow, title, accent)
    return s, y


def divider(kicker, title, sub):
    """Section divider. On the BMGF template the slate background, and the
    logo strip if the layout carries one, come from the layout itself; we drop
    its title placeholder and set our own type so the kicker/title/subtitle
    stack matches the rest of the deck."""
    if LAY_DIVIDER is None:
        s = prs.slides.add_slide(BLANK)
        rect(s, Inches(0), Inches(0), W, H, DARK)
    else:
        s = prs.slides.add_slide(LAY_DIVIDER)
        for ph in list(s.placeholders):
            ph._element.getparent().remove(ph._element)
    rect(s, MARGIN, Inches(2.85), Inches(1.5), Inches(0.05), DARK_KICKER)
    textbox(s, MARGIN, Inches(2.30), BODY_W, Inches(0.3),
            [(kicker.upper(), 13, True, DARK_KICKER)])
    textbox(s, MARGIN, Inches(3.15), BODY_W, Inches(1.0),
            [(title, 36, True, WHITE)], spacing=0.95, font=HEAD_FONT)
    textbox(s, MARGIN, Inches(4.35), Inches(9.2), Inches(1.0),
            [(sub, 15, False, DARK_SUB)], spacing=1.3)
    return s


def title_slide(main, sub, who, ver):
    """Opening slide. On the BMGF template this fills the branded title layout,
    so the parchment field, framing image and logo lockup are inherited."""
    if LAY_TITLE is None:
        s = prs.slides.add_slide(BLANK)
        rect(s, Inches(0), Inches(0), W, H, DARK)
        rect(s, MARGIN, Inches(2.28), Inches(2.0), Inches(0.06), DARK_KICKER)
        textbox(s, MARGIN, Inches(1.70), BODY_W, Inches(0.4),
                [("STARSIM / STISIM — HIVSIM ESWATINI", 14, True,
                  DARK_KICKER)])
        textbox(s, MARGIN, Inches(2.62), Inches(11.4), Inches(1.6),
                [(main, 44, True, WHITE), (sub, 21, False, DARK_SUB, 14)],
                spacing=0.98, font=HEAD_FONT)
        textbox(s, MARGIN, Inches(5.30), Inches(11.4), Inches(1.3),
                [(who, 15, True, RGBColor(0xE8, 0xE7, 0xE0)),
                 (ver, 13, False, RGBColor(0x9C, 0x9B, 0x92), 8)], spacing=1.3)
        return s
    # (text, colour, font). The date placeholder inherits white from the
    # layout, which is invisible on parchment, so every colour is set here.
    s = prs.slides.add_slide(LAY_TITLE)
    fills = {0: (main, INK, HEAD_FONT), 1: (sub, INK, BODY_FONT),
             14: (who, INK, BODY_FONT), 15: (ver, INK2, BODY_FONT)}
    for ph in list(s.placeholders):
        spec = fills.get(ph.placeholder_format.idx)
        if spec is None:
            ph._element.getparent().remove(ph._element)
            continue
        text, colour, face = spec
        ph.text_frame.text = text
        for para in ph.text_frame.paragraphs:
            for run in para.runs:
                run.font.name = face
                run.font.color.rgb = colour
    return s


def end_slide(title, sub):
    """Closing slide. Only exists on the BMGF template, which is where the
    logo lockup lives; returns None on the neutral theme."""
    if LAY_END is None:
        return None
    s = prs.slides.add_slide(LAY_END)
    textbox(s, MARGIN, Inches(3.05), Inches(10.0), Inches(1.2),
            [(title, 30, True, WHITE)], spacing=0.98, font=HEAD_FONT)
    textbox(s, MARGIN, Inches(4.15), Inches(9.4), Inches(1.0),
            [(sub, 14, False, DARK_SUB)], spacing=1.3)
    return s


def save(outfile):
    prs.save(outfile)
    # Output can live outside the repo (see build_story.py on why), so do not
    # assume relative_to(ROOT) succeeds.
    try:
        shown = outfile.relative_to(ROOT)
    except ValueError:
        shown = outfile
    print(f"wrote {shown}  ({len(prs.slides._sldIdLst)} slides)")
