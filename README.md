# Amyloid interface passivation model

This repository accompanies **Fibril interface passivation and clearance reshape bistability in a network model of amyloid associated damage**.

## Reviewed revision (15 October 2026)

- `paper_v3/Manuscript_Revised.pdf`: visually checked manuscript.
- `paper_v3/Manuscript_Revised.docx`: editable text; display equations are high-resolution images. Word pagination has not been visually validated because its renderer was unavailable.
- `paper_v3/Revision_and_Submission_Notes.pdf`: Chinese revision and journal memo.
- `manuscript/final_draft.md`: assembled source; `revised_body.md` is the reviewed authoring source.
- `manuscript/references.bib`, `reference_audit.json`, `references_verified.json`: bibliography and Crossref/PubMed verification records.
- `results_v3/data`: figure inputs, verified folds and matched-dose controls.
- `results_v3/figures`: vector PDF figures; Figs. 1, 4, 6 and S1 additionally include SVG.

The original `results_v2/data` are retained as inputs and for provenance. Legacy `profile_likelihood_ci.csv` is excluded from the revised manuscript: the old rerun script held nuisance parameters at reference values. The corrected rerun API now retains all six fitted parameters, but this revision makes no validated profile-likelihood or coverage claim.

## Model and interpretation

Each region carries monomer, oligomer and fibril burdens, free passivator concentration, interface occupancy and a reversible damage index. Free-interface activity suppresses both secondary nucleation and elongation. These are phenomenological burdens; secondary nucleation has no matching monomer sink, so total peptide mass is not stoichiometrically conserved. Row-normalized mixing preserves the strength-weighted sum on undirected raw graphs, not the arithmetic sum. Damage spreads from more damaged neighbours via the positive part of their damage difference.

For the reference q=1 model, analytic fast-variable elimination resolves the entire physical damage interval. The node folds are 0.928210 and 0.990198; the frozen-monomer threshold 0.547 is a reduced-system comparator, not an exact full-model treatment threshold. Two stable equilibria coexist between the folds. At occupancy 0.95, long integration reaches damage 0.0202 or 0.4499 depending on initial state.

The combination grid changes oligomer clearance and passivator infusion only; the treatment-window/Table 2 combination also enhances fibril clearance. Exact matched controls give mean Bliss excess 0.1529 and 90% positive points at tolerance 1e-6. Baseline-clearance points have zero excess within solver accuracy. Neither simulated synergy nor a low endpoint at time 200 establishes clinical efficacy or asymptotic network arrest.

## Reproduction

Install Python dependencies with `python -m pip install -r requirements.txt`. The recorded package versions are in `requirements-revision-lock.txt`. Simulation code works independently of the document tools. Current document builders use Cambria and SimSun fonts from Windows.

Fast rebuild from the supplied data:

```text
python reproduce_revision.py
python -m unittest discover -s tests -v
```

Recompute the affected numerical checks and exact single-agent controls as well:

```text
python reproduce_revision.py --recompute
```

This uses the archived original experiment data rather than rerunning the entire original fitting and 400-sample sensitivity ensemble. `run_pipeline.py` remains the original full-experiment driver. New reference checks can be refreshed with `python verify_references.py`; the build itself uses the local verified cache and needs no network access.

`publish_revision.py` defaults to a read-only upload plan. Only `--publish` updates the repository, as a verified fast-forward commit. It never uploads `.env` or credentials and preserves unselected remote paths. The existing nested checkout and its staged changes are not reset or committed by this workflow.

## Evidence limits

There is no clinical cohort or external kinetic validation. The connectome is a deterministic surrogate. Local Gaussian fitting intervals and the 12 independent marginal parameter draws are descriptive, not validated confidence or prediction intervals. Original results are preserved; revised statements distinguish approximation, numerical evidence and biological interpretation. Author details remain placeholders by request. No specific funding and no competing interests were reported.
