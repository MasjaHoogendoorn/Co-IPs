# Co-IP proteomics

Computational analysis of co-immunoprecipitation mass spectrometry (Co-IP-MS) data from *Plasmodium falciparum*.

Analyses by Masja Hoogendoorn, Swiss TPH.

This repository contains two independent experiments, each with its own analysis notebook, figures, results and interpretation.

## FR002 — PfSUN1 Co-IP

**Aim:** Identify proteins enriched in PfSUN1 pull-downs compared with controls.

| | |
|---|---|
| Bait | PfSUN1 (`Q8I5Q5`) |
| Samples | 3 pull-downs and 3 controls |
| Replication | 3 independent biological preparations per condition |
| Analysis | limma, proDA and SAINTexpress |

The analysis includes quality control, differential protein recovery and candidate interaction scoring. Biological replicates allow statistical comparisons, although the small sample size and missing measurements limit confidence.

**[FR002 analysis notebook](notebooks/FR002_SUN.ipynb) · [Results and interpretation](output/FR002/README.md)**

## CK001 — RON4 Co-IP

**Aim:** Compare proteins recovered with RON4 after DMSO and rapamycin treatment.

| | |
|---|---|
| Bait | RON4 (`Q8IIK5`) |
| Samples | DMSO, rapamycin and untagged NF54 control |
| Replication | 1 biological preparation per condition, measured twice |
| Analysis | Descriptive intensity comparisons and fold changes |

The analysis includes quality control, bait recovery and changes in co-recovered proteins. The repeated measurements are technical replicates, so differences are descriptive and no statistical significance testing is performed.

**[CK001 analysis notebook](notebooks/CK001_RON4.ipynb) · [Results and interpretation](output/CK001/README.md)**

## Installation

Create and activate the software environment:

```bash
mamba env create -f environment.yml -n coip
conda activate coip
```

The environment includes Python, R, limma and proDA.

## Input data

Both experiments use processed ProteoFlux tables, stored in `input/FR002/` and `input/CK001/`. The original workbooks are not modified.

**FR002 only:** SAINTexpress requires additional protein reference sequences. Download and prepare them using:

```bash
python input/reference/build_species_proteomes.py input/FR002/*.xlsx
```

FR002 also downloads SAINTexpress automatically on first use. Both steps require an internet connection.

CK001 requires neither step.

## Running an analysis

Start JupyterLab:

```bash
jupyter lab
```

Open the relevant notebook in `notebooks/` and run all cells from top to bottom.

The notebooks run independently, using shared Python and R functions from `src/coip/`. Experiment-specific settings are defined near the beginning of each notebook.

## Repository structure

```text
notebooks/       Analysis notebooks for FR002 and CK001
src/coip/        Shared Python and R code
input/           ProteoFlux tables and reference files
output/FR002/    PfSUN1 results
output/CK001/    RON4 results
environment.yml  Software dependencies
figure_style.txt Figure conventions
```

Each output directory contains figures (PDF, SVG and 400 dpi PNG), numerical result tables and an experiment-specific README. Re-running an analysis overwrites its generated outputs.

## Interpretation

Missing intensity measurements are kept as missing values, without imputation.

Co-IP-MS measures proteins recovered during a pull-down. Differences in recovery do not necessarily reflect changes in cellular protein abundance, and co-recovery does not establish direct physical binding.
