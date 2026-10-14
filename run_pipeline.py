"""Run the full experiment pipeline for the manuscript (v2)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ad_graph_model.experiments import run_all


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all v2 experiments.")
    parser.add_argument("--config", default="configs/v2_baseline.json")
    parser.add_argument("--outdir", default="results_v2")
    parser.add_argument("--skip-sensitivity", action="store_true")
    parser.add_argument("--n-sens", type=int, default=400)
    args = parser.parse_args()

    run_all(args.config, args.outdir,
            skip_sensitivity=args.skip_sensitivity, n_sens=args.n_sens)


if __name__ == "__main__":
    main()
