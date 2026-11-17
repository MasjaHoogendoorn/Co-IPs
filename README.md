# Co-IP proteomics

Co-immunoprecipitation mass spectrometry from *Plasmodium falciparum*, analysed by
Masja Hoogendoorn at Swiss TPH.

A Co-IP measures how much of a protein was recovered together with the bait. It does
not measure how much of that protein the parasite contains, and it is not evidence of
direct physical binding.

## The two experiments

| | FR002 | CK001 |
|---|---|---|
| bait | PfSUN1, `Q8I5Q5` | RON4, `Q8IIK5` |
| conditions | control and PfSUN1 pull-down | DMSO, rapamycin, untagged NF54 |
| runs per condition | 3 | 2 |
| replicate type | independent biological preparations | one preparation, injected twice |
| what that supports | limma and proDA, with p-values | descriptive fold changes only |
| specificity analysis | SAINTexpress on spectral counts | none |

CK001's two runs per condition are technical replicates, so it reports fold changes and
the measurements behind them and applies no statistical test. FR002 has three
independent preparations per condition, so it can.

`output/FR002/README.md` and `output/CK001/README.md` give the design, the main
findings and the limitations of each experiment. Read those before the notebooks.

## Software

```bash
mamba env create -f environment.yml -n coip
conda activate coip
```

That installs Python, R, limma, proDA and the rest. FR002 also uses SAINTexpress,
which the notebook downloads into `.tools/` on first use (about 9 MB) and checks
against its published checksum, so it needs an internet connection the first time.

## Input data

The two ProteoFlux workbooks are in the repository:

```
input/FR002/list_of_identified_proteins_FR002_proteoflux_table.xlsx
input/CK001/list_of_identified_proteins.xlsx
```

Each notebook records the SHA-256 of the workbook it read, and never modifies it.

FR002 also needs protein sequences, to give SAINTexpress the length of every prey. Two
of the three files are too large for the repository and are rebuilt from UniProt:

```bash
python input/reference/build_species_proteomes.py input/FR002/*.xlsx
```

This downloads the reviewed human proteome, the *P. falciparum* proteome and the
archived versions of accessions that UniProt has since retired. It takes a few minutes
and needs an internet connection. CK001 does not use these files.

## Running an analysis

Open either notebook in Jupyter and run all cells from top to bottom:

```bash
jupyter lab notebooks/FR002_SUN.ipynb
jupyter lab notebooks/CK001_RON4.ipynb
```

Each notebook is self-contained: it reads its workbook, writes its own figures and
tables, and depends on no other notebook. Paths are resolved relative to the project,
so the working directory does not matter.

Every scientific choice — the samples, the bait, the normalisation, the thresholds —
is in one configuration cell near the top. Nothing in `src/coip/` knows about a
particular experiment.

## Where the results go

```
output/FR002/figures/    11 figures as PDF, SVG and 400 dpi PNG
output/FR002/tables/     the source data behind every figure
output/FR002/saint/      the files SAINTexpress was given and what it returned
output/CK001/figures/    9 figures
output/CK001/tables/     the source data behind every figure
```

Re-running a notebook overwrites them.

## What the code does

```
src/coip/data.py          read and validate a workbook, QC numbers
src/coip/processing.py    log2, median equalisation, condition summaries
src/coip/differential.py  the Python side of limma and proDA
src/coip/limma.R          the limma model
src/coip/proda.R          the proDA model
src/coip/saint.py         spectral count summaries, and SAINTexpress for FR002
src/coip/figures.py       the figure style and every figure
```

`figure_style.txt` records the plotting conventions and the approved colours.

## Limitations that apply to both experiments

Missing measurements are never imputed and never drawn as zero. A protein that was not
quantified in a run counts as missing, and a fold change that cannot be estimated is
left empty rather than filled in.

A difference in recovery is not a difference in abundance, and agreement between an
intensity-based and a spectral-count analysis is not independent validation: both read
the same six measurements.
