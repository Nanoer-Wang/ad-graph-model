"""Convert the Stage-1 manuscript Markdown into a typeset HTML page.

This implements the doc-typeset `academic-paper` workflow:
  design tokens (academic-paper) -> :root CSS variables -> semantic
  <section role="cover"> + <section role="body"> -> figure/table markup.

The manuscript is an English-language SCI submission, so the GB/T 7713
Chinese numbering (一、二、三) is NOT applied: the manuscript keeps its own
"1. Introduction / 2. Model formulation / ..." numbering, and body text is
set left-aligned (ragged right) following English journal convention.
"""
from __future__ import annotations

import html
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "output" / "ad_passivation_paper" / "stage1" / "final_draft.md"
OUT_DIR = ROOT / "output" / "ad_passivation_paper" / "stage2"
OUT_HTML = OUT_DIR / "formatted-fibril-interface-passivation.html"


# --------------------------------------------------------------------------
# inline markdown -> HTML
# --------------------------------------------------------------------------
def inline(text: str) -> str:
    s = html.escape(text, quote=False)
    s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s, flags=re.S)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", s)
    return s


def split_row(line: str) -> list[str]:
    body = line.strip()
    if body.startswith("|"):
        body = body[1:]
    if body.endswith("|"):
        body = body[:-1]
    return [c.strip() for c in body.split("|")]


def is_separator(cells: list[str]) -> bool:
    return all(re.fullmatch(r":?-{2,}:?", c.replace(" ", "")) for c in cells if c)


# --------------------------------------------------------------------------
# block parsing
# --------------------------------------------------------------------------
def parse(md: str) -> tuple[str, list[str], list[tuple[str, str]]]:
    """Return (title, cover_paragraphs, body_blocks).

    body_blocks is a list of (kind, html) pairs.
    """
    lines = md.split("\n")
    title = ""
    cover: list[str] = []
    i = 0

    # --- title (first h1) -------------------------------------------------
    while i < len(lines) and not lines[i].startswith("# "):
        i += 1
    title = lines[i][2:].strip()
    i += 1

    # --- cover block: everything before the first h2 ----------------------
    cover: list[str] = []
    while i < len(lines) and not lines[i].startswith("## "):
        if lines[i].strip():
            cover.append(lines[i].strip())
        i += 1

    # --- body -------------------------------------------------------------
    blocks: list[tuple[str, str]] = []
    para: list[str] = []
    table_rows: list[list[str]] = []
    table_caption = ""
    ref_items: list[str] = []
    in_refs = False
    in_tables = False
    in_abstract = False

    def flush_para() -> None:
        nonlocal para
        if para:
            blocks.append(("p", " ".join(para)))
            para = []

    def flush_table() -> None:
        nonlocal table_rows, table_caption
        if not table_rows:
            return
        head, body_rows = table_rows[0], table_rows[1:]
        cap = (f"<caption>{inline(table_caption)}</caption>"
               if table_caption else "")
        th = "".join(f"<th>{inline(c)}</th>" for c in head)
        trs = "".join(
            "<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>"
            for r in body_rows)
        blocks.append(("table",
                       f'<table class="three-line-table">{cap}'
                       f"<thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>"))
        table_rows = []
        table_caption = ""

    def flush_refs() -> None:
        nonlocal ref_items, in_refs
        if ref_items:
            lis = "".join(f'<li class="reference-item">{inline(t)}</li>'
                          for t in ref_items)
            blocks.append(("refs", f'<ol class="reference-list">{lis}</ol>'))
            ref_items = []
        in_refs = False

    while i < len(lines):
        raw = lines[i].rstrip()
        s = raw.strip()

        if not s:
            if in_refs and ref_items:
                pass  # keep accumulating (multi-line refs unsupported)
            flush_para()
            i += 1
            continue

        if s.startswith("### "):
            flush_para()
            if in_refs:
                flush_refs()
            blocks.append(("h3", inline(s[4:].strip())))
            i += 1
            continue

        if s.startswith("## "):
            flush_para()
            flush_table()
            if in_refs:
                flush_refs()
            head = s[3:].strip()
            in_refs = head.lower().startswith("references")
            in_tables = head.lower().startswith("tables")
            in_abstract = head.lower().startswith("abstract")
            blocks.append(("h2", inline(head)))
            i += 1
            continue

        # images
        m = re.fullmatch(r"!\[(.*?)\]\((.*?)\)", s)
        if m:
            flush_para()
            alt, src = m.group(1), m.group(2)
            if alt.startswith("eq"):
                blocks.append(("eq", f'<img src="{src}" alt="{alt}">'))
            else:
                blocks.append(("figimg", f'<img src="{src}" alt="{alt}">'))
            i += 1
            continue

        # ordered-list reference items
        if in_refs and re.match(r"^\d+\.\s", s):
            ref_items.append(re.sub(r"^\d+\.\s*", "", s))
            i += 1
            continue

        # figure / table captions
        m = re.match(r"^\*\*(Figure|Table)\s+(\d+)\.\*\*\s*(.*)$", s)
        if m:
            flush_para()
            flush_table()
            kind, num, rest = m.group(1), m.group(2), m.group(3)
            cap = f"{kind} {num}. {rest}"
            if kind == "Figure":
                blocks.append(("figcap", inline(cap)))
            else:
                table_caption = cap
            i += 1
            continue

        # table rows
        if s.startswith("|"):
            flush_para()
            cells = split_row(s)
            if not is_separator(cells):
                table_rows.append(cells)
            i += 1
            continue

        # keywords line
        if s.startswith("**Keywords:**"):
            flush_para()
            blocks.append(("keywords", inline(s)))
            i += 1
            continue

        # ordinary paragraph
        para.append(s)
        i += 1

    flush_para()
    flush_table()
    flush_refs()
    return title, cover, blocks


# --------------------------------------------------------------------------
# HTML assembly
# --------------------------------------------------------------------------
CSS = """
  @page { @bottom-center { content: counter(page); } }
  @page cover { @bottom-center { content: none; } }
  section[role="cover"] { page: cover; }

  :root {
    --ff-heading: var(--typography-fontFamily-bodyLatin, "Times New Roman", serif);
    --ff-body: var(--typography-fontFamily-bodyLatin, "Times New Roman", serif);
    --fs-cover-title: var(--typography-fontSize-title, 18pt);
    --fs-cover-meta: var(--typography-fontSize-h3, 12pt);
    --fs-h2: var(--typography-fontSize-h2, 14pt);
    --fs-h3: var(--typography-fontSize-h3, 12pt);
    --fs-h4: var(--typography-fontSize-body, 12pt);
    --fs-body: var(--typography-fontSize-body, 12pt);
    --fs-small: var(--typography-fontSize-caption, 10.5pt);
    --lh-body: var(--typography-lineHeight-body, 1.5);
    --lh-heading: var(--typography-lineHeight-heading, 1.5);
    --fw-bold: var(--typography-fontWeight-heading, bold);
    --fw-normal: var(--typography-fontWeight-body, normal);
    --color-primary: var(--color-primary, #000000);
    --color-text: var(--color-text, #000000);
    --color-muted: var(--color-text, #000000);
    --color-border: var(--color-tableBorder, #000000);
    --color-bg: var(--color-background, #ffffff);
    --spacing-paragraph: var(--spacing-paragraph, 0.5em);
    --spacing-section: var(--spacing-sectionGap, 1.5em);
    --spacing-block: var(--spacing-paragraph, 0.5em);
    --margin-page: var(--layout-marginTop, 2.5cm);
    --page-content-width: 16.0cm;
  }

  body {
    font-family: var(--ff-body);
    font-size: var(--fs-body);
    line-height: var(--lh-body);
    color: var(--color-text);
    background: var(--color-bg);
    max-width: var(--page-content-width);
    margin-left: auto;
    margin-right: auto;
  }

  h2 {
    font-family: var(--ff-heading); font-size: var(--fs-h2);
    font-weight: var(--fw-bold); line-height: var(--lh-heading);
    margin-top: var(--spacing-section); margin-bottom: var(--spacing-paragraph);
    text-align: left; color: var(--color-text);
  }
  h3 {
    font-family: var(--ff-heading); font-size: var(--fs-h3);
    font-weight: var(--fw-bold); line-height: var(--lh-heading);
    margin-top: var(--spacing-block); margin-bottom: var(--spacing-paragraph);
    text-align: left; color: var(--color-text);
  }
  h4 {
    font-family: var(--ff-heading); font-size: var(--fs-h4);
    font-weight: var(--fw-bold);
    margin-top: var(--spacing-block); margin-bottom: var(--spacing-paragraph);
    text-align: left; color: var(--color-text);
  }
  p {
    margin-bottom: 0; text-align: justify;
    font-family: var(--ff-body); font-size: var(--fs-body);
    color: var(--color-text); text-indent: 2em; line-height: 1.45;
  }
  /* journal convention: the first paragraph after a heading, figure,
     equation or table is set flush left (no first-line indent) */
  h2 + p, h3 + p, p.caption + p, p.figure-img + p, p.equation + p,
  .three-line-table + p, .abstract-text p, .abstract-keywords p {
    text-indent: 0; margin-bottom: var(--spacing-paragraph);
  }

  .cover-line {
    text-align: center; font-size: var(--fs-cover-meta);
    margin-bottom: var(--spacing-paragraph); color: var(--color-text);
  }
  .cover-note {
    font-size: var(--fs-small); margin-top: var(--spacing-section);
    text-align: center; color: var(--color-text);
  }
  .cover-eyebrow {
    text-align: center; font-size: var(--fs-small);
    margin-top: var(--spacing-section); color: var(--color-text);
  }
  .cover-title {
    text-align: center; font-size: var(--fs-cover-title);
    font-weight: var(--fw-bold); font-family: var(--ff-heading);
    margin-top: var(--spacing-section); margin-bottom: var(--spacing-section);
    color: var(--color-text); line-height: var(--lh-heading);
  }
  .cover-line {
    text-align: center; font-size: var(--fs-cover-meta);
    margin-bottom: var(--spacing-paragraph); color: var(--color-text);
  }

  /* abstract block: border-left => requires leading &nbsp;&nbsp; */
  .abstract {
    margin-bottom: var(--spacing-section); padding: var(--spacing-block);
    border-left: 3px solid var(--color-border);
  }
  .abstract-heading {
    text-align: center; font-weight: var(--fw-bold); font-size: var(--fs-body);
    margin-bottom: var(--spacing-paragraph); color: var(--color-text);
  }
  .abstract-text { text-align: left; font-size: var(--fs-body); }
  .abstract-keywords { text-align: left; font-size: var(--fs-body); }

  /* figures and equations */
  p.figure-img { text-align: center; margin: var(--spacing-block) 0; }
  p.figure-img img { max-width: 100%; }
  p.equation { text-align: center; margin: var(--spacing-block) 0; }
  p.equation img { max-width: 100%; }
  p.caption {
    font-size: var(--fs-small); text-align: justify;
    margin-top: 0.4em; margin-bottom: var(--spacing-section);
    color: var(--color-text); text-indent: 0;
  }

  /* three-line table */
  .three-line-table {
    border-collapse: collapse; width: 100%;
    margin-top: var(--spacing-block); margin-bottom: var(--spacing-section);
  }
  .three-line-table caption {
    font-size: var(--fs-small); margin-bottom: 0.4em; text-align: left;
  }
  .three-line-table thead tr:first-child { border-top: 2px solid var(--color-border); }
  .three-line-table thead tr:last-child { border-bottom: 1px solid var(--color-border); }
  .three-line-table tbody tr:last-child { border-bottom: 2px solid var(--color-border); }
  .three-line-table td, .three-line-table th {
    border-left: none; border-right: none;
    padding: 0.3em 0.6em; text-align: left; font-size: var(--fs-small);
  }
  .three-line-table th { font-weight: var(--fw-bold); background: transparent; }

  /* references */
  .references-block { margin-top: var(--spacing-section); }
  .reference-list { padding-left: 2em; }
  .reference-item {
    margin-bottom: 0.4em; font-size: var(--fs-small); text-align: left;
  }
"""


def build(title: str, cover: list[str], blocks) -> str:
    out: list[str] = []
    out.append('<!DOCTYPE html>')
    out.append('<html lang="en">')
    out.append('<head>')
    out.append('<meta charset="UTF-8">')
    out.append('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
    out.append('<meta name="docx-page-size" content="A4">')
    out.append(f"<title>{html.escape(title)}</title>")
    out.append("<style>")
    out.append(CSS)
    out.append("</style>")
    out.append("</head>")
    out.append("<body>")

    def clean_cover(line: str) -> str:
        """Strip the markdown italic/bold markers used in the placeholder block."""
        s = line.strip()
        s = re.sub(r"^\*\*(.*?)\*\*", r"\1", s)
        s = re.sub(r"^\*Placeholder format:\*\s*", "", s)
        s = s.strip("*").strip()
        s = re.sub(r"^\*\s*(Correspondence:)", r"\1", s)
        return s.strip().rstrip(",").strip()

    def cover_class(line: str) -> str:
        if line.startswith("**Running head:**") or "Running head" in line:
            return "cover-eyebrow"
        if "Author A. Name" in line:
            return "cover-line"
        if line.lstrip().startswith(("\u00b9", "\u00b2", "\u00b3")):
            return "cover-line"
        if "Correspondence" in line:
            return "cover-line"
        return "cover-line cover-note"

    # ---------------- cover ----------------
    out.append('  <section role="cover">')
    out.append(f'    <h1 class="cover-title">{html.escape(title)}</h1>')
    ordered = ([l for l in cover if "Author A. Name" in l]
               + [l for l in cover
                  if l.lstrip().startswith(("\u00b9", "\u00b2", "\u00b3"))]
               + [l for l in cover if "Correspondence" in l]
               + [l for l in cover if cover_class(l) == "cover-line cover-note"]
               + [l for l in cover if cover_class(l) == "cover-eyebrow"])
    for line in ordered:
        out.append(f'    <p class="{cover_class(line)}">'
                   f'{html.escape(clean_cover(line))}</p>')
    out.append("  </section>")

    # ---------------- body ----------------
    out.append('  <section role="body" data-page-restart="1">')

    in_abstract = False
    pending_fig: str | None = None
    for kind, payload in blocks:
        if kind == "h2":
            head = payload
            if re.match(r"^Abstract", head, re.I):
                out.append('    <div class="abstract" aria-label="Abstract">')
                out.append('      <p class="abstract-heading">Abstract</p>')
                in_abstract = True
                continue
            if in_abstract:
                out.append("    </div>")
                in_abstract = False
            out.append(f"    <h2>{head}</h2>")
            continue

        if kind == "h3":
            out.append(f"    <h3>{payload}</h3>")
            continue

        if kind == "figimg":
            pending_fig = payload
            continue

        if kind == "figcap":
            if pending_fig:
                out.append(f'    <p class="figure-img">{pending_fig}</p>')
                pending_fig = None
            out.append(f'    <p class="caption">{payload}</p>')
            continue

        if kind == "eq":
            out.append(f'    <p class="equation">{payload}</p>')
            continue

        if kind == "keywords":
            out.append(f'      <p class="abstract-keywords">{payload}</p>')
            continue

        if kind == "table":
            out.append(f"    {payload}")
            continue

        if kind == "refs":
            out.append(f'    <div class="references-block" aria-label="References">')
            out.append(f"      {payload}")
            out.append("    </div>")
            continue

        # ordinary paragraph
        if in_abstract:
            out.append(f'      <p class="abstract-text">&nbsp;&nbsp;{payload}</p>')
        else:
            out.append(f"    <p>{payload}</p>")

    if in_abstract:
        out.append("    </div>")
    out.append("  </section>")
    out.append("</body>")
    out.append("</html>")
    return "\n".join(out) + "\n"


FIG_WIDTH_CM = 16.0      # A4 text width at 2.5 cm margins
EQ_MAX_WIDTH_CM = 15.0   # leave a little air around long equations


def annotate_images(html_text: str, img_dir: Path) -> str:
    """Add explicit width + centring style to every <img>.

    html4docx opens ``src`` relative to the process CWD and, when no width
    is given, inserts the picture at its native (DPI-derived) size — the
    600 dpi figures would then be ~26 cm wide on a 16 cm text block.  Widths
    are therefore pinned here.
    """
    cache: dict[str, float] = {}

    def natural_cm(name: str) -> float:
        if name not in cache:
            from PIL import Image
            with Image.open(img_dir / name) as im:
                dpi = im.info.get("dpi", (96, 96))
                dx = dpi[0] if isinstance(dpi, (tuple, list)) else dpi
                cache[name] = im.size[0] / float(dx) * 2.54
        return cache[name]

    def repl(m: re.Match) -> str:
        src = m.group("src")
        name = src.split("/")[-1]
        nat = natural_cm(name)
        if name.startswith("eq"):
            w = min(nat, EQ_MAX_WIDTH_CM)
        else:
            w = min(nat, FIG_WIDTH_CM)
        return (f'<img src="{src}" width="{w:.2f}cm" '
                f'style="display:block;margin-left:auto;margin-right:auto;">')

    return re.sub(r'<img src="(?P<src>[^"]+)"[^>]*>', repl, html_text)


def main() -> None:
    md = SRC.read_text(encoding="utf-8")
    title, cover, blocks = parse(md)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "intermediate").mkdir(parents=True, exist_ok=True)

    # images used by the HTML must sit next to it
    img_src = SRC.parent / "images"
    img_dst = OUT_DIR / "images"
    if img_dst.exists():
        shutil.rmtree(img_dst)
    shutil.copytree(img_src, img_dst)

    html_text = annotate_images(build(title, cover, blocks), img_dst)
    OUT_HTML.write_text(html_text, encoding="utf-8")

    counts: dict[str, int] = {}
    for kind, _ in blocks:
        counts[kind] = counts.get(kind, 0) + 1
    print("title:", title[:70], "...")
    print("cover lines:", len(cover))
    print("blocks:", counts)
    print("images copied:", len(list(img_dst.glob("*.png"))))
    print("html:", OUT_HTML)


if __name__ == "__main__":
    main()
