"""Embed the eight manuscript figures inline (figure + caption) at the
position of their first substantive discussion, and drop the standalone
captions section.  Idempotent: safe to re-run (checks for existing embeds).
"""
from __future__ import annotations

import re
from pathlib import Path

DRAFT = Path(__file__).resolve().parent / \
    "output" / "ad_passivation_paper" / "stage1" / "final_draft.md"

# figure number -> (anchor: unique tail of the paragraph after which the
# figure is inserted, citation to add inside that paragraph or None)
TARGETS = {
    1: ("which we use throughout as the model's direct measure of the "
        "passivation effect.", None),
    2: ("a self-reinforcing loop that terminates only when neuronal loss "
        "drives production to zero at a \u2192 1.", None),
    3: ("an important robustness result for a model whose connectome is a "
        "surrogate.", None),
    4: ("the intermediate band is the feedback gap, in which the subsystem "
        "is stable but the full system is not.", None),
    5: ("achieves at the reference dose what monotherapy requires a large "
        "overdose to achieve (Section 4.5).", None),
    6: ("the quantitative basis for combination regimens in which each "
        "agent is used at a dose it could not have been used at alone.",
        None),
    7: ("affect neither timing nor outcome at the network level.", None),
    8: ("which is the prerequisite for future calibration against "
        "biomarker data.", None),
}


def main() -> None:
    text = DRAFT.read_text(encoding="utf-8")

    # ---- extract the standalone captions
    cap_pat = re.compile(r"\*\*Figure (\d)\.\*\* (.+?)(?=\n\n|\Z)", re.S)
    caps = {int(m.group(1)): m.group(2).strip()
            for m in cap_pat.finditer(text)}
    missing = set(range(1, 9)) - set(caps)
    if missing:
        raise SystemExit(f"captions missing: {sorted(missing)}")

    # ---- remove the standalone captions section (last section in file)
    i = text.find("## Figure captions")
    if i < 0:
        raise SystemExit("captions section not found")
    if text.find("## ", i + 5) >= 0:
        raise SystemExit("captions section is not last; aborting")
    text = text[:i].rstrip("\n") + "\n"
    print("captions section removed")

    # ---- insert figures inline
    for n in range(1, 9):
        if f"![fig{n}](" in text:
            print(f"fig{n} already embedded, skipping")
            continue
        anchor = TARGETS[n][0]
        idx = text.find(anchor)
        if idx < 0:
            raise SystemExit(f"anchor for fig{n} not found")
        end = idx + len(anchor)
        block = (f"\n\n![fig{n}](images/fig{n}.png)\n\n"
                 f"**Figure {n}.** {caps[n]}")
        text = text[:end] + block + text[end:]
        print(f"fig{n} embedded")

    # ---- cite Fig. 1 in the model-formulation text
    if "(Fig. 1)" not in text:
        old = ("the network- and time-averaged aggregation efficiency is")
        new = ("the network- and time-averaged aggregation efficiency is")
        # add citation at the end of the 'Two features' paragraph instead
        anchor = ("the transition from healthy to damaged behaviour is "
                  "graded at the single-node level and sharp only through "
                  "the nonlinear feedback of the network.")
        i = text.find(anchor)
        if i >= 0:
            text = text[:i + len(anchor)] + \
                " Fig. 1 summarises the model architecture." + text[i + len(anchor):]
            print("Fig. 1 citation added")

    DRAFT.write_text(text, encoding="utf-8")
    print("done")


if __name__ == "__main__":
    main()
