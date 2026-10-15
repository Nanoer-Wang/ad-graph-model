"""Validate and package the reviewed article and reproducible source files.

Only explicitly selected files enter the archive. Local credentials, checkouts,
temporary files and unrelated original document drafts are never selected.
"""
from pathlib import Path
import hashlib
import json
import re
import zipfile

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'output' / 'Submission_Package_900dpi.zip'
PREFIX = 'Submission_Package_900dpi/'
FIGURES = [
    'fig1_model_schematic', 'fig2_baseline_dynamics',
    'fig3_spatial_propagation', 'fig4_bifurcation_bistability',
    'fig5_passivation_dose_response', 'fig6_therapy_design',
    'fig7_global_sensitivity', 'fig8_identifiability', 'figS1_fold_detail',
]
SCRIPTS = [
    'make_figures.py', 'run_pipeline.py', 'revision_analysis.py',
    'revision_bliss.py', 'verify_references.py', 'assemble_revision.py',
    'render_equations.py', 'md_to_html_paper.py', 'build_pdf.py',
    'build_manuscript_docx.py', 'build_docx.py', 'make_revision_notes.py',
    'rerun_profile_ci.py', 'analysis_network_stats.py', 'reproduce_revision.py',
    'package_submission.py', 'publish_revision.py',
]
README = '''投稿材料包

Manuscript：英文正文的 DOCX 与 PDF。正文两端对齐；文中引用及文末条目采用 [1] 数字顺序格式，文末使用悬挂缩进。
Figures_900dpi：8张正文图和1张补充图，均从绘图代码原生渲染为900 dpi，提供PNG和无损LZW压缩RGB TIFF；并非提高旧图的DPI标签或插值放大。
Figures_vector：同一组图的矢量PDF，以及可用的SVG。排版需要缩放时优先采用矢量文件。
References：37篇引用的BibTeX、核对记录和元数据缓存。
Notes：修订说明、科学证据边界及期刊建议。
Code：模型、绘图与文章生成代码，配置、原始归档和修订结果数据、源稿、依赖版本及回归检查。README.md说明复现方法。
MANIFEST.json：每个文件的字节数与SHA-256；同时记录图片像素、分辨率、稿件页数及排版检查结果。

使用方法
1. 投稿系统的正文文件选用 Manuscript 内的 DOCX 或 PDF。
2. 独立图像上传使用 Figures_900dpi/TIFF，图号对应 Figure_1 至 Figure_8 及 Figure_S1。
3. 解压后进入 Code，安装 requirements.txt 中的依赖，运行 python reproduce_revision.py 即可从归档数据重建。完整原始实验入口为 run_pipeline.py。文章生成需要Windows的Cambria、SimSun字体；模型分析不受此限制。
4. 作者、单位和通讯信息按用户要求保留占位，基金与利益冲突已声明为无。正式递交前补齐作者信息，并按目标期刊指南调整专属格式。不同SCI期刊的引用体例并不统一，本包采用数字顺序Vancouver体例。

核验范围
PDF已逐页渲染检查；DOCX已检查正文对齐、37个参考文献标签、公式与图片嵌入、表格结构。当前环境缺少可用的Word渲染器，先前用于准备渲染器的操作被自动审批拒绝，因此DOCX的实际分页仍需在Word中复核。DOCX与PDF采用同一正文源，二者分页不保证完全相同。展示公式以900 dpi图片嵌入；正文可编辑，公式原始表达式随代码提供。
本包不包含.env、GitHub token、Git历史或本地临时文件。
'''


def main():
    files = {}

    def add(path, destination):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(path)
        if destination in files:
            raise ValueError('Duplicate archive destination: ' + destination)
        files[destination] = path

    docx = ROOT / 'paper_v3/Manuscript_Revised.docx'
    pdf = ROOT / 'paper_v3/Manuscript_Revised.pdf'
    doc = Document(docx)
    assert doc.styles['Normal'].paragraph_format.alignment == WD_ALIGN_PARAGRAPH.JUSTIFY
    refs = [p for p in doc.paragraphs if p.paragraph_format.left_indent
            and p.paragraph_format.first_line_indent
            and p.paragraph_format.first_line_indent < 0]
    assert len(refs) == 37
    for n, para in enumerate(refs, 1):
        assert para.text.startswith(f'[{n}] '), para.text[:40]
    reader = PdfReader(pdf)
    pdftext = '\n'.join(p.extract_text() for p in reader.pages)
    assert all(f'[{n}]' in pdftext for n in range(1, 38))
    assert len(doc.inline_shapes) == 24
    add(docx, 'Manuscript/Manuscript.docx')
    add(pdf, 'Manuscript/Manuscript.pdf')

    image_records = []
    for i, stem in enumerate(FIGURES, 1):
        label = f'Figure_{i}' if i <= 8 else 'Figure_S1'
        figdir = ROOT / 'results_v3/figures'
        for ext, folder in [('png', 'PNG'), ('tif', 'TIFF')]:
            path = figdir / (stem + '.' + ext)
            with Image.open(path) as im:
                # Pillow defaults to a warning above 89 MP; these native
                # scientific plots are known, trusted locally generated files.
                dpi = im.info.get('dpi', (0, 0))
                assert min(dpi) >= 899.9, (path, dpi)
                assert min(im.size) > 3000, (path, im.size)
                if ext == 'tif':
                    assert im.mode == 'RGB'
                    assert im.info.get('compression') == 'tiff_lzw'
                image_records.append(dict(figure=label, format=ext,
                    pixels=list(im.size), dpi=[float(v) for v in dpi], color_mode=im.mode))
                im.verify()
            add(path, f'Figures_900dpi/{folder}/{label}.{ext}')
        for ext in ['pdf', 'svg']:
            path = figdir / (stem + '.' + ext)
            if path.exists():
                add(path, f'Figures_vector/{label}.{ext}')

    for name in ['references.bib', 'reference_audit.json', 'references_verified.json']:
        add(ROOT / 'manuscript' / name, 'References/' + name)
    for name in ['Revision_and_Submission_Notes.pdf', 'REVISION_NOTES.txt']:
        add(ROOT / 'paper_v3' / name, 'Notes/' + name)
    for name in SCRIPTS + ['README.md', 'requirements.txt', 'requirements-revision-lock.txt']:
        add(ROOT / name, 'Code/' + name)
    for pattern in ['src/ad_graph_model/*.py', 'configs/*.json', 'tests/*.py',
                    'manuscript/*.md', 'manuscript/*.bib', 'manuscript/*.json',
                    'manuscript/images/*.png', 'results_v2/data/*', 'results_v3/data/*']:
        for path in sorted(ROOT.glob(pattern)):
            if path.is_file():
                add(path, 'Code/' + path.relative_to(ROOT).as_posix())
    source = 'output/ad_passivation_paper/stage1/final_draft.md'
    add(ROOT / source, 'Code/' + source)

    # Credential bytes are checked in memory, with no logging of their value.
    token = None
    env = ROOT / '.env'
    if env.exists():
        for line in env.read_text(encoding='utf-8').splitlines():
            if '=' in line:
                key, val = line.split('=', 1)
                if key.strip() == 'GIT_WANG':
                    token = val.strip().strip('\"').strip("'").encode()
    manifest = dict(manuscript_pdf_pages=len(reader.pages), references=37,
        docx_normal_alignment='justified', reference_labels='[1] ... [37]',
        docx_visual_render_verified=False, images=image_records, files=[])
    OUT.parent.mkdir(exist_ok=True)
    temp = OUT.with_suffix('.zip.partial')
    with zipfile.ZipFile(temp, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for rel, path in sorted(files.items()):
            assert not any(s in rel.split('/') for s in ['.git', '.env', '__pycache__'])
            raw = path.read_bytes()
            if token and token in raw:
                raise RuntimeError('Credential bytes found in a selected file')
            manifest['files'].append(dict(path=rel, bytes=len(raw),
                sha256=hashlib.sha256(raw).hexdigest()))
            z.writestr(PREFIX + rel, raw)
        raw = README.encode('utf-8-sig')
        manifest['files'].append(dict(path='README_投稿包说明.txt', bytes=len(raw),
            sha256=hashlib.sha256(raw).hexdigest()))
        z.writestr(PREFIX + 'README_投稿包说明.txt', raw)
        z.writestr(PREFIX + 'MANIFEST.json', json.dumps(manifest, ensure_ascii=False, indent=2))
    with zipfile.ZipFile(temp) as z:
        assert z.testzip() is None
        for record in manifest['files']:
            assert hashlib.sha256(z.read(PREFIX + record['path'])).hexdigest() == record['sha256']
    temp.replace(OUT)
    print(json.dumps(dict(archive=str(OUT), bytes=OUT.stat().st_size,
        files=len(manifest['files']) + 1, figures=len(FIGURES),
        pdf_pages=len(reader.pages)), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
