# Delayed recall and subsequent ADL disability in HAALSI and CLHLS

Python analysis and figure code accompanying the manuscript prepared for **BMC Geriatrics**.
Release date: 19 October 2026. This folder is the repository root: upload its **contents** to GitHub.

The primary outcome contains five shared activities: bathing, dressing, toileting, transfers and eating. Analyses are performed separately for HAALSI and CLHLS. No pooled effect or meta-analysis is calculated.

## What this release contains

- `audit/independent_audit.py`: self-contained raw-file reconstruction, complete-case primary and sensitivity models, cohort flow counts, and HC0 numerical checks. This is the recommended main analysis entry point, reached via `run_analysis.py main`.
- `analysis/five_item_MI.py`: original statsmodels FCS implementation, seed **20261011**, 40 imputations and 20 cycles. It is retained unchanged from the archived original analysis, and requires the extra dependency file.
- `audit/score_quality_and_MI.py`: HAALSI raw-word-validity sensitivities and independent MI checks, seed **20261017**.
- `audit/rerun_original_HA_MI.py`: independent HAALSI released-score MI check, seed **20261018**.
- `scripts/plot_forest.py`: Figure 1, reading the full-precision aggregate results rather than rounded manuscript entries. It can instead read freshly reconstructed model JSON using `--results`.
- `reference_results/audited_results_full_precision.csv`: archived aggregate estimates with explicit version labels, including both markers. These are reference results, not newly estimated values.
- `analysis/*.json`, `audit/*results.json`, `flow_pvalues/*.json`: aggregate reference/check outputs. There are no participant-level records in this release. The main audit compares a fresh raw-data reconstruction against these saved reference results; this is intentional numerical regression checking, not estimation from saved effects.

## Install

Use Python **3.12** (the archived independent audit used **3.12.14**).

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell, use instead:
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

To rerun the original statsmodels MI implementation as well:

```bash
python -m pip install -r requirements-original-mi.txt
```

NumPy, pandas and SciPy are pinned to the independently audited environment. Matplotlib is pinned to the environment used for Figure 1 in this release. The exact version of statsmodels used by the original seed-20261011 run was not retained in the recovered audit metadata; the extra dependency file therefore specifies a compatible range rather than claiming an exact historical lock. Keep the installed version in your own execution log. The independent NumPy/SciPy raw-data reconstruction does not require statsmodels.

## Checks and plotting without raw data

```bash
python run_analysis.py verify-archive
python run_analysis.py figure
```

The first command checks 38 archived effect/P-value rows, checks all four primary figure estimates against the independent audit, and tests the solver against an analytic intercept-only logistic result and the missing/positive outcome rule. These checks do **not** rerun raw-data models. The second command reproduces the forest plot from archived aggregate estimates and writes PNG, TIFF, PDF, SVG and source CSV under `figures/`.

## Obtain the raw data yourself

The original individual-level data are **not included** and should not be added to the public repository. Acquire the authorised releases from their providers and comply with registration, access and redistribution terms.

| Local filename under `inputs/` | Provider | Required archive member |
|---|---|---|
| `HAALSI W1-3 Longitudinal Dataset.dta` | [HAALSI W1–3, Harvard Dataverse](https://doi.org/10.7910/DVN/Q2JFPV) | Direct `.dta` file |
| `dataverse_files.zip` | [HAALSI Wave 4, Harvard Dataverse](https://doi.org/10.7910/DVN/JNQCPB) | `HAALSI Wave 4 Public Dataset.dta` |
| `clhls_2018_cross_sectional_dataset_15874.zip` | [CLHLS 2018, PKU](https://doi.org/10.18170/DVN/WBO7LK) | The 2018 cross-sectional `.dta`, not the 1998–2018 combined longitudinal file |
| `dataverse_files (1).zip` | [CLHLS 2021, PKU](https://doi.org/10.18170/DVN/UNYEL1) | `CLHLS2021_released.dta` |

These are the input filenames and formats used by the archived code. If a current provider release supplies `.tab` rather than `.dta`, it is not interchangeable by renaming the extension: obtain the matching Stata release or implement and validate a typed reader before running this code. ZIP members must have the names shown above at archive root. The 2018 ZIP reader selects its first `.dta`; keep only the intended 2018 cross-sectional data file in that archive.

`inputs/input_manifest.json` records the four original SHA-256 hashes and sizes. Downloaded ZIP containers may have different hashes even when the dataset is unchanged; inspect release versions and members when this occurs. By default a mismatch warns; `--strict-hashes` makes it an error.

## Run the raw-data analysis

```bash
python run_analysis.py check-inputs
python run_analysis.py main
# Optional original MI, seed 20261011; reconstructs data before running MI
python run_analysis.py mi-original
# Optional independent score-validity / fresh-seed MI checks
python run_analysis.py mi-audit
# All stages; requires the statsmodels extra dependency
python run_analysis.py all
```

Each raw-data command checks that the four inputs exist. `main` reconstructs variables from raw files, fits the complete-case models, compares them with archived references and generates the forest plot from **fresh** model results. The main output is `audit/independent_audit_results.json`. Model failures or reference discrepancies raise errors rather than silently using archived estimates. Small numerical differences and fresh-seed MI differences must be evaluated according to the tolerances and methods, not selected by P value.

Optional alternative data location:

```bash
# Linux/macOS
export HAALSI_CLHLS_INPUTS=/absolute/path/to/authorised/input-files
# Windows PowerShell
# $env:HAALSI_CLHLS_INPUTS = 'C:\path\to\authorised\input-files'
python run_analysis.py main
```

BLAS threads are set to one for child analysis processes. Raw-data stages create temporary **participant-level** checkpoints, `audit/private_rebuilt_frames.pkl` and (for original MI) `analysis/five_item_frames.pkl`, on the user's computer. They are excluded by `.gitignore`; do not publish them. New aggregate outputs may overwrite archived aggregate JSON files in a local clone: retain the tagged release as the fixed comparison version.

## Expected primary results

| Cohort | n / outcomes | Poorer recall OR (95% CI) | Psychological symptoms OR (95% CI) |
|---|---:|---|---|
| HAALSI | 1199 / 105 | 1.328 (1.065–1.655) | 1.064 (0.878–1.291) |
| CLHLS | 3824 / 623 | 1.173 (1.078–1.277) | 1.114 (1.009–1.229) |

Both markers enter the same cohort-specific model. Units are one cohort-specific baseline reference SD, calculated among baseline-five-activity-negative participants with both scores observed, before outcome/covariate complete-case restrictions. Complete-case inference uses HC0 sandwich SEs and two-sided nominal Wald z tests. The age term is `(age − 75) / 10` plus its square. Other adjustment variables are sex, any schooling, hypertension, diabetes, cardiovascular/cerebrovascular disease, smoking, drinking and partner status. MI pools HC0 completed-data variances using Rubin/Barnard–Rubin inference.

## Measurement and interpretation

- HAALSI captures help received; CLHLS captures help requirements/dependence. The shared domains do not establish complete measurement equivalence.
- Any valid positive activity defines a positive composite; a negative composite requires all included activities to be valid negatives; other patterns remain unknown.
- HAALSI eating screen-negative/missing-help responses are treated as structural nonreceipt. Wave 4 documentation supports this; exact matching Wave 2 routing documentation remains a limitation.
- The primary released HAALSI recall score is retained. The conservative raw-word-validity sensitivity excludes two original-main records, yielding **1197 / 104**. That restriction is not applied to the primary figure.
- The sixth activity differs: indoor walking in HAALSI, continence in CLHLS. Six-activity analyses are sensitivities, not the primary outcome.
- Death, unknown follow-up and uncompleted cognitive scores are not imputed. Models are conditional follow-up-status associations, not competing-risk models or causal estimates. All P values are nominal.

## Validation and provenance of this release

The 17 October archive contains evidence of raw-data reconstruction and numerical agreement of complete-case models, and independent MI checks. On **19 October**, this release was assembled from those archived scripts and aggregate results. The scripts were syntax-checked, archived estimates checked, the analytic solver/outcome-logic checks passed and Figure 1 regenerated. **A new end-to-end raw-data rerun was not performed during this packaging step because the four raw inputs were not present in the current workspace.** The original statsmodels MI code was recovered and syntax-checked; it was not rerun in this packaging environment.

The numerical core of the independent audit is retained. Its only code adaptation is support for `HAALSI_CLHLS_INPUTS`; new files provide orchestration, checks, documentation and plotting. The original seed-20261011 MI implementation is unchanged. This is a code/results release; the Word manuscript is supplied separately.

## Uploading to GitHub

1. Create your repository and upload the contents of this directory (including `.gitignore`).
2. Keep raw data and participant-level checkpoints out of the repository.
3. Review the README, contribution statements and an appropriate code license with the authors. No license has been selected on their behalf.
4. Tag the code version used for submission. A Zenodo archive can provide a permanent DOI for the tagged release.
5. Insert the **actual** public repository URL/DOI in the manuscript after publication of the repository. This ZIP has not itself been uploaded to GitHub or Zenodo.

Before the code is shared as the final reproducibility release, an author with the four authorised inputs should run the raw-data stages in the intended environment and retain the execution logs. This does not imply a change to the manuscript's existing estimates.
