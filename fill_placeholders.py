"""Fill the remaining manuscript placeholders (Table 4 and Section 4.7)
from the recomputed fit / profile-CI outputs.

Run AFTER rerun_profile_ci.py has refreshed profile_likelihood_ci.csv.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "results_v2" / "data"
DRAFT = ROOT / "output" / "ad_passivation_paper" / "stage1" / "final_draft.md"

DISPLAY = {"kn": "kn", "ks": "ks", "ke": "ke", "cY": "cY",
           "etaY": "\u03b7Y", "KY": "KY"}


def fmt(v: float, digits: int = 3) -> str:
    """Significant-figure formatting that keeps small values readable."""
    if v == 0:
        return "0"
    from math import log10, floor
    exp = floor(log10(abs(v)))
    dec = max(digits - 1 - exp, 0)
    if exp >= digits:
        return f"{v:.0f}"
    return f"{v:.{dec}f}"


def main() -> None:
    fit = json.load(open(DATA / "fit_result.json"))
    ci = pd.read_csv(DATA / "profile_likelihood_ci.csv").set_index("parameter")
    rec = fit["recovery"]

    text = DRAFT.read_text(encoding="utf-8")

    # ---- Table 4 rows
    for i, name in enumerate(["kn", "ks", "ke", "cY", "etaY", "KY"], start=1):
        r = ci.loc[name]
        true_v = rec[name]["true"]
        best_v = rec[name]["fit"]
        rel_err = rec[name]["rel_err"]
        row = (f"| {DISPLAY[name]} | {fmt(true_v)} | {fmt(best_v)} | "
               f"[{fmt(r['ci_lo'])}, {fmt(r['ci_hi'])}] | "
               f"{rel_err * 100:+.1f}% |")
        pat = re.compile(rf"\| \[\[ID{i}\]\] \| \[\[ID{i}_T\]\] \| \[\[ID{i}_B\]\] "
                         rf"\| \[\[ID{i}_CI\]\] \| \[\[ID{i}_E\]\] \|")
        text, n = pat.subn(row, text)
        if n != 1:
            raise SystemExit(f"row {i} not found (n={n})")

    # ---- Section 4.7 scalars
    max_err = fit["max_rel_err"] * 100.0
    text = text.replace("[[IDENT_MAXERR]]", f"{max_err:.1f}")

    def rel_width(name: str) -> float:
        r = ci.loc[name]
        return float((r["ci_hi"] - r["ci_lo"]) / r["best"])

    text = text.replace("[[CI_KS]]", f"{rel_width('ks') * 100:.0f}%")
    text = text.replace("[[CI_CY]]", f"{rel_width('cY') * 100:.0f}%")

    leftover = re.findall(r"\[\[[A-Z0-9_]+\]\]", text)
    if leftover:
        print("WARNING leftover placeholders:", sorted(set(leftover)))
    else:
        print("all placeholders filled")

    DRAFT.write_text(text, encoding="utf-8")
    print(ci.to_string())
    print(f"max recovery error: {max_err:.2f}%")
    print(f"relative CI widths: ks {rel_width('ks') * 100:.1f}%, "
          f"cY {rel_width('cY') * 100:.1f}%")


if __name__ == "__main__":
    main()
