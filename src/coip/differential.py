"""Run R limma on a log2 intensity matrix.

limma is the reference implementation of the moderated t-test, so the
statistics are done by limma itself rather than reimplemented here. This
module writes the matrix and the design to temporary CSV files, calls
src/coip/limma.R with Rscript, and reads the result back.

Nothing about a particular experiment lives here: the design is an argument.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

LIMMA_SCRIPT = Path(__file__).with_name("limma.R")


def find_rscript(explicit=None):
    """Locate an Rscript that has limma installed.

    Looks at the argument, then $COIP_RSCRIPT, then PATH, then the other conda
    environments next to the running interpreter , which is how R usually
    reaches a Python kernel on a shared machine. Raises with instructions if
    none of them works.
    """
    candidates = [explicit, os.environ.get("COIP_RSCRIPT"), shutil.which("Rscript")]
    candidates += sorted(str(p) for p in
                         Path(sys.prefix).parent.glob("*/bin/Rscript"))
    for candidate in candidates:
        if candidate and Path(candidate).exists() and has_limma(candidate):
            return str(candidate)
    raise RuntimeError(
        "No Rscript with limma was found. Install it with\n"
        "  conda install -c conda-forge r-base bioconductor-limma\n"
        "or, in Colab,\n"
        "  !apt-get -qq install r-base && "
        "Rscript -e 'install.packages(\"BiocManager\"); "
        "BiocManager::install(\"limma\")'\n"
        "then set COIP_RSCRIPT to its path.")


def has_limma(rscript) -> bool:
    try:
        done = subprocess.run(
            [str(rscript), "-e",
             'cat(requireNamespace("limma", quietly=TRUE))'],
            capture_output=True, text=True, timeout=180)
    except (OSError, subprocess.SubprocessError):
        return False
    return done.returncode == 0 and "TRUE" in done.stdout


def limma_contrast(log2_intensity, groups, treatment, reference, rscript=None,
                   keep_files=None):
    """Moderated t-test of `treatment` against `reference`.

    log2_intensity  protein groups in rows, samples in columns, no missing
                    values (filter to complete cases first)
    groups          sample -> group label, for every column used
    rscript         path to an Rscript with limma; discovered if omitted
    keep_files      directory to leave the limma input/output CSVs in, for
                    inspection; a temporary directory is used otherwise

    Returns a DataFrame indexed like the input, with log2FC, avg_log2, t,
    p_value, p_adjusted, B and the degrees of freedom the test used.
    """
    rscript = rscript or find_rscript()

    samples = [s for s in log2_intensity.columns if str(s) in groups]
    unknown = [s for s in log2_intensity.columns if str(s) not in groups]
    if unknown:
        raise KeyError(f"no group given for sample(s) {unknown}")
    matrix = log2_intensity[samples]
    if matrix.isna().any().any():
        raise ValueError("the matrix still has missing values; "
                         "filter to complete cases before testing")
    for level in (treatment, reference):
        if sum(groups[str(s)] == level for s in samples) < 2:
            raise ValueError(f"group '{level}' has fewer than 2 samples; "
                             "limma cannot estimate a residual variance")

    workdir = Path(keep_files) if keep_files else Path(tempfile.mkdtemp(prefix="coip_limma_"))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        intensity_file = workdir / "limma_input_log2_intensities.csv"
        design_file = workdir / "limma_input_design.csv"
        out_file = workdir / "limma_output.csv"

        named = matrix.copy()
        named.insert(0, "protein_id", [f"row{i}" for i in matrix.index])
        named.to_csv(intensity_file, index=False)
        pd.DataFrame({"sample": [str(s) for s in samples],
                      "group": [groups[str(s)] for s in samples]}
                     ).to_csv(design_file, index=False)

        done = subprocess.run(
            [str(rscript), str(LIMMA_SCRIPT), str(intensity_file), str(design_file),
             str(treatment), str(reference), str(out_file)],
            capture_output=True, text=True)
        if done.returncode != 0:
            raise RuntimeError(f"limma.R failed:\n{done.stdout}\n{done.stderr}")

        result = pd.read_csv(out_file)
        result.index = matrix.index
        result = result.drop(columns=["protein_id"])
        result.attrs["limma_message"] = done.stdout.strip()
        return result
    finally:
        if not keep_files:
            shutil.rmtree(workdir, ignore_errors=True)


PRODA_SCRIPT = Path(__file__).with_name("proda.R")


def proda_contrast(log2_intensity, groups, treatment, reference, rscript=None,
                   keep_files=None):
    """Differential abundance with proDA, keeping the missing values.

    Where `limma_contrast` needs complete cases, this takes the matrix as it
    is: a protein not measured in a sample stays NA and contributes through
    proDA's dropout model. Nothing is imputed, so a protein seen only in the
    pull-downs can still be estimated , with the uncertainty that implies,
    which is reported as `se`, `df` and `n_obs` rather than hidden.

    proDA does not normalise between samples, so the result is on whatever
    scale the input was on.

    Returns a DataFrame indexed like the input.
    """
    rscript = rscript or find_rscript()

    samples = [s for s in log2_intensity.columns if str(s) in groups]
    unknown = [s for s in log2_intensity.columns if str(s) not in groups]
    if unknown:
        raise KeyError(f"no group given for sample(s) {unknown}")
    matrix = log2_intensity[samples]
    for level in (treatment, reference):
        if sum(groups[str(s)] == level for s in samples) < 2:
            raise ValueError(f"group '{level}' has fewer than 2 samples")
    measured = matrix.notna().sum(axis=1)
    if (measured == 0).any():
        raise ValueError(f"{int((measured == 0).sum())} proteins have no "
                         "measurement at all; proDA cannot place them")

    workdir = Path(keep_files) if keep_files else Path(
        tempfile.mkdtemp(prefix="coip_proda_"))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        intensity_file = workdir / "proda_input_log2_intensities.csv"
        design_file = workdir / "proda_input_design.csv"
        out_file = workdir / "proda_output.csv"

        named = matrix.copy()
        named.insert(0, "protein_id", [f"row{i}" for i in matrix.index])
        named.to_csv(intensity_file, index=False)
        pd.DataFrame({"sample": [str(s) for s in samples],
                      "group": [groups[str(s)] for s in samples]}
                     ).to_csv(design_file, index=False)

        done = subprocess.run(
            [str(rscript), str(PRODA_SCRIPT), str(intensity_file),
             str(design_file), str(treatment), str(reference), str(out_file)],
            capture_output=True, text=True)
        if done.returncode != 0:
            raise RuntimeError(f"proda.R failed:\n{done.stdout}\n{done.stderr}")

        result = pd.read_csv(out_file)
        result.index = matrix.index
        result = result.drop(columns=["protein_id"])
        result.attrs["proda_message"] = done.stdout.strip()
        result.attrs["proda_warnings"] = done.stderr.strip()
        return result
    finally:
        if not keep_files:
            shutil.rmtree(workdir, ignore_errors=True)
