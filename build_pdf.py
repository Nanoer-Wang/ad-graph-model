"""Export the manuscript as a submission-quality PDF.

No Word / LibreOffice is installed on this machine and Edge's headless
``--print-to-pdf`` is unreliable here (it silently produces nothing whenever
a browser instance is already running), so the PDF is laid out directly with
reportlab from the parsed Stage-1 Markdown.

Typography notes:
  * Cambria is used because it is the only widely available Windows serif on
    this machine that covers every glyph used in the manuscript (Greek
    letters, ∈, ≫, ≲, ≳, →, ⟨ ⟩) — Times New Roman misses seven of them.
  * Page numbers are stamped bottom-centre; the title page is unnumbered and
    the body restarts at 1, matching the .docx.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from md_to_html_paper import SRC, inline, parse  # noqa: E402

from reportlab.lib import colors  # noqa: E402
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle  # noqa: E402
from reportlab.lib.units import cm  # noqa: E402
from reportlab.pdfbase import pdfmetrics  # noqa: E402
from reportlab.pdfbase.ttfonts import TTFont  # noqa: E402
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,  # noqa: E402
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

PDF = ROOT / "paper_v3" / "Manuscript_Revised.pdf"
TEXT_WIDTH = 16.0 * cm          # A4 minus 2.5 cm margins
MARGIN = 2.5 * cm

FONTS = Path("C:/Windows/Fonts")
BODY = "Cambria"


# --------------------------------------------------------------------------
# fonts
# --------------------------------------------------------------------------
def register_fonts() -> None:
    regular = FONTS / "cambria.ttc"
    pdfmetrics.registerFont(TTFont(BODY, str(regular), subfontIndex=0))
    pdfmetrics.registerFont(TTFont(BODY + "-Bold", str(FONTS / "cambriab.ttf")))
    pdfmetrics.registerFont(TTFont(BODY + "-Italic", str(FONTS / "cambriai.ttf")))
    pdfmetrics.registerFont(TTFont(BODY + "-BoldItalic",
                                   str(FONTS / "cambriaz.ttf")))
    from reportlab.lib.fonts import addMapping
    for style, name in ((0, BODY), (1, BODY + "-Bold"),
                        (2, BODY + "-Italic"), (3, BODY + "-BoldItalic")):
        addMapping(BODY, style, 0, name)


# --------------------------------------------------------------------------
# styles
# --------------------------------------------------------------------------
def styles() -> dict[str, ParagraphStyle]:
    s = {}
    s["title"] = ParagraphStyle("title", fontName=BODY + "-Bold", fontSize=16,
                                leading=21, alignment=TA_CENTER,
                                spaceAfter=6)
    s["cover"] = ParagraphStyle("cover", fontName=BODY, fontSize=10.5,
                                leading=15, alignment=TA_CENTER)
    s["cover-note"] = ParagraphStyle("cover-note", fontName=BODY + "-Italic",
                                     fontSize=9, leading=13,
                                     alignment=TA_CENTER, spaceBefore=10)
    s["h2"] = ParagraphStyle("h2", fontName=BODY + "-Bold", fontSize=12.5,
                             leading=16, alignment=TA_LEFT,
                             spaceBefore=16, spaceAfter=7)
    s["h3"] = ParagraphStyle("h3", fontName=BODY + "-Bold", fontSize=11,
                             leading=14, alignment=TA_LEFT,
                             spaceBefore=12, spaceAfter=5)
    s['h2'].keepWithNext = True
    s['h3'].keepWithNext = True
    s["body"] = ParagraphStyle("body", fontName=BODY, fontSize=11,
                               leading=15.5, alignment=TA_JUSTIFY,
                               spaceAfter=0, firstLineIndent=0)
    # Journal convention: consecutive body paragraphs are separated by a
    # first-line indent rather than by vertical space.  The first paragraph
    # after a heading, figure, table or equation is set flush left (style
    # "body"); later ones use "body-indent".
    s["body-indent"] = ParagraphStyle("body-indent", parent=s["body"],
                                      firstLineIndent=1.1 * 10.5)
    s["abstract"] = ParagraphStyle("abstract", fontName=BODY, fontSize=10,
                                   leading=14.5, alignment=TA_JUSTIFY,
                                   leftIndent=1.1 * cm, rightIndent=1.1 * cm,
                                   spaceAfter=6)
    s["abstract-head"] = ParagraphStyle("abstract-head",
                                        fontName=BODY + "-Bold", fontSize=11,
                                        alignment=TA_CENTER, spaceAfter=8,
                                        spaceBefore=4)
    s["keywords"] = ParagraphStyle("keywords", fontName=BODY, fontSize=10,
                                   leading=14, alignment=TA_JUSTIFY,
                                   leftIndent=1.1 * cm, rightIndent=1.1 * cm,
                                   spaceBefore=4, spaceAfter=10)
    s["caption"] = ParagraphStyle("caption", fontName=BODY, fontSize=9,
                                  leading=12, alignment=TA_JUSTIFY,
                                  spaceBefore=4, spaceAfter=12)
    s["tablecap"] = ParagraphStyle("tablecap", fontName=BODY, fontSize=9,
                                   leading=12, alignment=TA_LEFT,
                                   spaceBefore=10, spaceAfter=3)
    s["cell"] = ParagraphStyle("cell", fontName=BODY, fontSize=8.5, leading=11)
    s["cellhead"] = ParagraphStyle("cellhead", fontName=BODY + "-Bold",
                                   fontSize=8.5, leading=11)
    s["ref"] = ParagraphStyle("ref", fontName=BODY, fontSize=9, leading=12.5,
                              alignment=TA_LEFT, leftIndent=0.8 * cm,
                              firstLineIndent=-0.8 * cm, spaceAfter=3)
    return s


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def esc(text: str) -> str:
    return (text.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;"))


def clean_cover(line: str) -> str:
    """Strip the markdown markers used in the placeholder author block."""
    s = line.strip()
    s = re.sub(r"^\*\*(.*?)\*\*", r"\1", s)
    s = re.sub(r"^\*Placeholder format:\*\s*", "", s)
    s = s.strip("*").strip()
    s = re.sub(r"^\*\s*(Correspondence:)", r"\1", s)
    return s.strip().rstrip(",").strip()


def ordered_cover(cover: list[str]) -> list[str]:
    authors = [c for c in cover if "Author A. Name" in c]
    affils = [c for c in cover
              if c.lstrip().startswith(("\u00b9", "\u00b2", "\u00b3"))]
    corr = [c for c in cover if "Correspondence" in c]
    notes = [c for c in cover
             if c not in authors + affils + corr and "Running head" not in c]
    heads = [c for c in cover if "Running head" in c]
    return authors + affils + corr + notes + heads


def rl(text: str) -> str:
    """HTML produced by md_to_html_paper -> reportlab paragraph markup."""
    t = (text.replace("<strong>", "<b>").replace("</strong>", "</b>")
             .replace("<em>", "<i>").replace("</em>", "</i>"))
    return t


def natural_cm(path: Path) -> float:
    from PIL import Image as PILImage
    with PILImage.open(path) as im:
        dpi = im.info.get("dpi", (96, 96))
        dx = dpi[0] if isinstance(dpi, (tuple, list)) else dpi
        return im.size[0] / float(dx) * 2.54


IMG_DIR = SRC.parent / "images"


def image_size(name: str) -> tuple[float, float]:
    """Display size (cm) for an image, honouring its native aspect ratio.

    Both dimensions are returned explicitly: reportlab otherwise derives the
    height itself and, for high-resolution figures, can end up with a box taller
    than the text block (LayoutError).
    """
    from PIL import Image as PILImage
    with PILImage.open(IMG_DIR / name) as im:
        w_px, h_px = im.size
        dpi = im.info.get("dpi", (96, 96))
        dx = dpi[0] if isinstance(dpi, (tuple, list)) else dpi
        dy = dpi[1] if isinstance(dpi, (tuple, list)) else dx
        w_cm, h_cm = w_px / float(dx) * 2.54, h_px / float(dy) * 2.54
    max_w = 14.8 if name.startswith("eq") else 15.9
    max_h = 19.2
    scale = min(1.0, max_w / w_cm, max_h / h_cm)
    return w_cm * scale, h_cm * scale


def rows_from_html(html_text: str) -> tuple[str, list[list[str]]]:
    m = re.search(r"<caption>(.*?)</caption>", html_text, re.S)
    caption = m.group(1) if m else ""
    body = re.sub(r"<caption>.*?</caption>", "", html_text, flags=re.S)
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", body, re.S):
        cells = re.findall(r"<t[hd]>(.*?)</t[hd]>", tr, re.S)
        rows.append(cells)
    return caption, rows


def col_widths(rows: list[list[str]], avail: float) -> list[float]:
    n = max(len(r) for r in rows)
    if rows[0][:2] == ['Symbol','Meaning']:
        return [avail*w for w in [.12,.31,.14,.43]]
    if rows[0][0] == 'Strategy':
        return [avail*w for w in [.36,.12,.16,.18,.18]]
    if rows[0][0] == 'Rank':
        return [avail*w for w in [.07,.50,.18,.25]]
    if rows[0][0] == 'Parameter' and n == 5:
        return [avail*w for w in [.16,.23,.19,.19,.23]]
    weight = [max((len(r[i]) if i < len(r) else 0) for r in rows) + 2
              for i in range(n)]
    total = sum(weight)
    return [w / total * avail for w in weight]


def make_table(rows: list[list[str]], st: dict[str, ParagraphStyle],
               avail: float) -> Table:
    widths = col_widths(rows, avail)
    data = []
    for ri, row in enumerate(rows):
        style = st["cellhead"] if ri == 0 else st["cell"]
        data.append([Paragraph(rl(c), style) for c in row])
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("LINEABOVE", (0, 0), (-1, 0), 1.4, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 1.4, colors.black),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]))
    return t


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------
def build_flow(title, cover, blocks, st):
    flow: list = []
    avail = TEXT_WIDTH

    # ---- cover ----
    flow.append(Spacer(1, 2.0 * cm))
    flow.append(Paragraph(rl(title), st["title"]))
    flow.append(Spacer(1, 0.6 * cm))
    for line in ordered_cover(cover):
        style = st["cover-note"] if ("to be completed" in line
                                     or "Running head" in line) else st["cover"]
        flow.append(Paragraph(esc(clean_cover(line)), style))
    flow.append(PageBreak())

    in_abstract = False
    equation_number = 0
    main_equations = []
    pending_img: str | None = None
    for i, (kind, payload) in enumerate(blocks):
        if kind == "h2":
            if payload == 'Supplementary material':
                flow.append(PageBreak())
            if in_abstract:
                in_abstract = False
            if re.match(r"^Abstract", payload, re.I):
                flow.append(Paragraph("Abstract", st["abstract-head"]))
                in_abstract = True
            else:
                flow.append(Paragraph(rl(payload), st["h2"]))
            continue

        if kind == "h3":
            flow.append(Paragraph(rl(payload), st["h3"]))
            continue

        if kind == "figimg":
            m = re.search(r'src="images/([^"]+)"', payload)
            pending_img = m.group(1) if m else None
            continue

        if kind == "figcap":
            items = []
            if pending_img:
                w, h = image_size(pending_img)
                items.append(Image(str(IMG_DIR / pending_img),
                                   width=w * cm, height=h * cm,
                                   hAlign="CENTER"))
                pending_img = None
            items.append(Paragraph(rl(payload), st["caption"]))
            flow.append(KeepTogether(items))
            continue

        if kind == "eq":
            m = re.search(r'src="images/([^"]+)"', payload)
            if m:
                w, h = image_size(m.group(1))
                equation_number += 1
                equation = Image(str(IMG_DIR / m.group(1)),width=w*cm,height=h*cm,hAlign="CENTER")
                row=Table([[equation,Paragraph(f'({equation_number})',st['body'])]],colWidths=[15.0*cm,1*cm])
                row.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'MIDDLE'),('ALIGN',(0,0),(0,0),'CENTER'),('LEFTPADDING',(0,0),(-1,-1),0),('RIGHTPADDING',(0,0),(-1,-1),0),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),8)]))
                if equation_number <= 6:
                    main_equations.append(row)
                    if equation_number == 6:
                        flow.append(KeepTogether(main_equations))
                else:
                    flow.append(row)
            continue

        if kind == "keywords":
            flow.append(Paragraph(rl(payload), st["keywords"]))
            continue

        if kind == "table":
            caption, rows = rows_from_html(payload)
            if not rows:
                continue
            parts = []
            if caption:
                parts.append(Paragraph(rl(caption), st["tablecap"]))
            parts.append(make_table(rows, st, avail))
            parts.append(Spacer(1, 6))
            flow.append(KeepTogether(parts))
            continue

        if kind == "refs":
            items = re.findall(r'<li class="reference-item">(.*?)</li>',
                               payload, re.S)
            for i, it in enumerate(items, 1):
                flow.append(Paragraph(f"[{i}] {rl(it)}", st["ref"]))
            continue

        # "p" blocks still carry raw markdown emphasis; every other kind was
        # already run through inline() by the parser
        prev_kind = blocks[i - 1][0] if i > 0 else None
        if not in_abstract and prev_kind not in ("h2", "h3", "figcap",
                                                 "table", "eq", "figimg"):
            style = st["body-indent"]
        else:
            style = st["abstract"] if in_abstract else st["body"]
        flow.append(Paragraph(rl(inline(payload)), style))

    return flow


def page_number(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont(BODY, 10)
    canvas.drawCentredString(A4[0] / 2.0, 1.25 * cm,
                             str(canvas.getPageNumber() - 1))
    canvas.restoreState()


def main() -> None:
    register_fonts()
    md = SRC.read_text(encoding="utf-8")
    title, cover, blocks = parse(md)
    st = styles()
    flow = build_flow(title, cover, blocks, st)

    PDF.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(str(PDF), pagesize=A4,
                          leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=MARGIN, bottomMargin=MARGIN,
                          title=title, author="Author A. Name; Author B. Name",
                          subject="Submitted manuscript")
    frame = Frame(MARGIN, MARGIN, TEXT_WIDTH, A4[1] - 2 * MARGIN, id="body",
                  leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates([
        PageTemplate(id="cover", frames=[frame]),
        PageTemplate(id="numbered", frames=[frame], onPage=page_number),
    ])
    # page 1 uses the first template (unnumbered title page); the explicit
    # PageBreak in the flow then switches every later page to "numbered"
    flow.insert(0, NextPageTemplate("numbered"))

    doc.build(flow)
    print(f"pdf:   {PDF} ({PDF.stat().st_size/1e6:.1f} MB)")
    print("pages:", _page_count(PDF))


def _page_count(path: Path) -> int:
    from pypdf import PdfReader
    return len(PdfReader(str(path)).pages)


if __name__ == "__main__":
    main()
