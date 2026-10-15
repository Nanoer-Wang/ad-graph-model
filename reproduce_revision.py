"""Rebuild the reviewed manuscript, optionally recomputing revised analyses."""
from pathlib import Path
import argparse,subprocess,sys
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--recompute',action='store_true');args=ap.parse_args()
    scripts=(['revision_analysis.py','revision_bliss.py'] if args.recompute else [])+[
        'make_figures.py','assemble_revision.py','render_equations.py','build_pdf.py','build_manuscript_docx.py','make_revision_notes.py']
    for script in scripts:
        extra=['--outdir','manuscript/images'] if script=='render_equations.py' else []
        subprocess.run([sys.executable,str(ROOT/script),*extra],cwd=ROOT,check=True)
if __name__=='__main__':main()
