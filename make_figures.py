"""Render all manuscript figures from the experiment data."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ad_graph_model.figures_core import (fig1_schematic, fig2_baseline,
                                         fig3_spatial, fig4_bifurcation)
from ad_graph_model.figures_therapy import (fig5_dose, fig6_therapy,
                                            fig7_sensitivity,
                                            fig8_identifiability)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render all figures.")
    parser.add_argument("--data", default="results_v2/data")
    parser.add_argument("--outdir", default="results_v2/figures")
    parser.add_argument("--only", default=None,
                        help="comma-separated figure numbers, e.g. 1,3,5")
    args = parser.parse_args()

    data = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    only = None
    if args.only:
        only = {int(x) for x in args.only.split(",")}

    jobs = {
        1: lambda: fig1_schematic(outdir, data / "baseline_timeseries.csv",
                                  data / "strategy_Passivation_timeseries.csv"),
        2: lambda: fig2_baseline(data, outdir),
        3: lambda: fig3_spatial(data, outdir),
        4: lambda: fig4_bifurcation(data, outdir),
        5: lambda: fig5_dose(data, outdir),
        6: lambda: fig6_therapy(data, outdir),
        7: lambda: fig7_sensitivity(data, outdir),
        8: lambda: fig8_identifiability(data, outdir),
    }
    for num in sorted(jobs):
        if only is not None and num not in only:
            continue
        paths = jobs[num]()
        print(f"Fig {num}: " + ", ".join(str(p) for p in paths), flush=True)


if __name__ == "__main__":
    main()
