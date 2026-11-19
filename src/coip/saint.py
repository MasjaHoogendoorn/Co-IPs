"""Score interaction confidence with SAINTexpress.

SAINTexpress compares the spectral counts a prey reaches in a bait pull-down
with the counts it reaches in the negative controls, and reports a probability
that the interaction is specific. The model is published (Teo et al. 2014,
J Proteomics 100:37-43) and we run the authors' own compiled program rather
than reimplementing it.

The program is a command-line tool: it reads three tab-separated files and
writes 'list.txt'. This module builds those three tables from a count matrix,
runs the program and reads the result back. Nothing here knows about a
particular experiment.
"""

from __future__ import annotations

import hashlib
import platform
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

# The release published on the saint-apms SourceForge project. The checksum is
# verified after download, so a tampered or truncated archive cannot be used.
SAINTEXPRESS_URL = ("https://sourceforge.net/projects/saint-apms/files/"
                    "SAINTexpress_v3.6.3__2018-03-09.tar.gz/download")
SAINTEXPRESS_SHA256 = \
    "668dada2f3617bc9bb899314c097bd41a9de8f9bed78e934459cc125126b1482"
SAINTEXPRESS_RELEASE = "SAINTexpress_v3.6.3__2018-03-09"

_PLATFORM_DIR = {("Linux", "x86_64"): "Linux64",
                 ("Windows", "AMD64"): "Windows64"}


def find_saintexpress(cache_dir, program="SAINTexpress-spc", explicit=None,
                      download=True):
    """Path to a SAINTexpress executable, downloading the release if needed.

    Looks at `explicit`, then at anything already unpacked under `cache_dir`,
    then at PATH, and only then fetches the published archive. Returns the
    path; raises if the platform has no precompiled binary.
    """
    cache_dir = Path(cache_dir)
    unpacked = cache_dir / SAINTEXPRESS_RELEASE
    key = (platform.system(), platform.machine())
    if key not in _PLATFORM_DIR:
        raise RuntimeError(
            f"SAINTexpress ships no precompiled binary for {key}. Build it from "
            f"the sources in {SAINTEXPRESS_URL} and pass the path explicitly.")
    binary = unpacked / "Precompiled_binaries" / _PLATFORM_DIR[key] / program

    for candidate in (explicit, binary if binary.exists() else None,
                      shutil.which(program)):
        if candidate and Path(candidate).exists():
            return Path(candidate)
    if not download:
        raise FileNotFoundError(f"{program} not found and download=False")

    cache_dir.mkdir(parents=True, exist_ok=True)
    archive = cache_dir / f"{SAINTEXPRESS_RELEASE}.tar.gz"
    if not archive.exists():
        urllib.request.urlretrieve(SAINTEXPRESS_URL, archive)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != SAINTEXPRESS_SHA256:
        archive.unlink()
        raise RuntimeError(f"checksum of the SAINTexpress archive is {digest}, "
                           f"expected {SAINTEXPRESS_SHA256}")
    with tarfile.open(archive) as tar:
        tar.extractall(cache_dir, filter="data")
    binary.chmod(0o755)
    for sibling in binary.parent.iterdir():
        sibling.chmod(0o755)
    return binary


def version(binary):
    done = subprocess.run([str(binary), "-v"], capture_output=True, text=True)
    return (done.stdout + done.stderr).strip()


# --- preys -------------------------------------------------------------------

def prey_table(protein_groups, gene_names=None):
    """One SAINTexpress prey row per protein group.

    A protein group holds several accessions that the search could not tell
    apart. SAINTexpress wants a single identifier, so the first accession
    becomes the prey id, and the group it stands for is reported alongside it
    so nothing is lost.

    No prey length is recorded. The prey file SAINTexpress-spc reads takes two
    columns, identifier and gene name, as well as three with a length between
    them, and the length never reaches the score: the program stores it,
    `get_preyLength` is never called, and 'list.txt' has no length column.
    Length normalisation belonged to the original SAINT, not to this one, so
    no protein sequences are needed to score a pull-down.

    Returns the prey table and a mapping back to the original protein groups.
    """
    rows = []
    for position, group in enumerate(protein_groups):
        members = [m.strip() for m in str(group).split(";") if m.strip()]
        rows.append({
            "row": position,
            "prey_id": members[0],
            "protein_group": str(group),
            "n_accessions": len(members),
            "gene_name": (str(gene_names[position]).split(";")[0]
                          if gene_names is not None else members[0]) or members[0],
        })
    mapping = pd.DataFrame(rows)
    if mapping["prey_id"].duplicated().any():
        duplicated = mapping.loc[mapping["prey_id"].duplicated(), "prey_id"].tolist()
        raise ValueError(f"prey ids are not unique: {duplicated[:5]}")
    return mapping


# --- input tables ------------------------------------------------------------

def _is_number(text):
    """Whether SAINTexpress would read `text` as a number rather than a name."""
    try:
        float(text)
    except ValueError:
        return False
    return True


def saint_tables(counts, prey_mapping, bait_of, control_samples):
    """Build SAINTexpress's three input tables from a spectral count matrix.

    counts           protein groups in rows, samples in columns, raw counts
    prey_mapping     as returned by prey_table
    bait_of          sample -> bait name (the pull-down it belongs to)
    control_samples  the samples that are negative controls

    Counts are passed through unchanged. SAINTexpress does its own
    normalisation internally and must see the original numbers.

    The prey file is written with two columns, identifier and gene name; see
    `prey_table` for why no length is needed.
    """
    samples = list(counts.columns)
    unknown = [s for s in samples if s not in bait_of]
    if unknown:
        raise KeyError(f"no bait assigned to sample(s) {unknown}")
    if not set(control_samples) <= set(samples):
        raise KeyError("control samples are not all in the count matrix")

    # SAINTexpress reads a second column it can parse as a number as a prey
    # length, and then takes the identifier as the gene name. A missing name
    # arrives here as 'nan', which parses. Either way the name would be lost
    # silently, so neither is written.
    names = prey_mapping["gene_name"].astype(str)
    unusable = [n for n in names if not n.strip() or _is_number(n)]
    if unusable:
        raise ValueError(
            f"{len(unusable)} prey gene name(s) are missing or would be read "
            f"as a prey length: {unusable[:5]}")
    prey = pd.DataFrame({"prey_id": prey_mapping["prey_id"],
                         "gene_name": names})

    bait = pd.DataFrame({
        "ip_id": samples,
        "bait_id": [bait_of[s] for s in samples],
        "control": ["C" if s in set(control_samples) else "T" for s in samples]})

    long = (counts.loc[prey_mapping["row"]]
            .set_axis(prey_mapping["prey_id"].to_numpy())
            .rename_axis("prey_id").reset_index()
            .melt(id_vars="prey_id", var_name="ip_id", value_name="count"))
    long = long[long["count"].notna() & (long["count"] > 0)]
    inter = pd.DataFrame({
        "ip_id": long["ip_id"],
        "bait_id": [bait_of[s] for s in long["ip_id"]],
        "prey_id": long["prey_id"],
        "count": long["count"].astype(int)}).sort_values(["ip_id", "prey_id"])

    return {"inter": inter.reset_index(drop=True), "prey": prey, "bait": bait}


def run(tables, binary, workdir=None, options=()):
    """Run SAINTexpress on prepared tables and return its result.

    `options` are passed to the program unchanged, e.g. ('-L2',) to cap the
    number of controls used per prey. The three input files and list.txt are
    left in `workdir` when one is given, so a run can be inspected or repeated
    by hand.
    """
    for name in ("inter", "prey", "bait"):
        if name not in tables:
            raise KeyError(f"missing the '{name}' table")

    temporary = workdir is None
    workdir = Path(workdir or tempfile.mkdtemp(prefix="coip_saint_"))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        for name in ("inter", "prey", "bait"):
            tables[name].to_csv(workdir / f"{name}.dat", sep="\t",
                                header=False, index=False)
        command = [str(Path(binary).resolve()), *options,
                   "inter.dat", "prey.dat", "bait.dat"]
        done = subprocess.run(command, cwd=workdir, capture_output=True, text=True)
        result_file = workdir / "list.txt"
        if done.returncode != 0 or not result_file.exists():
            raise RuntimeError(
                f"SAINTexpress failed ({done.returncode}):\n{done.stdout}\n{done.stderr}")
        result = pd.read_csv(result_file, sep="\t")
        result.attrs["command"] = " ".join(command)
        result.attrs["log"] = done.stdout.strip()
        return result
    finally:
        if temporary:
            shutil.rmtree(workdir, ignore_errors=True)


def reference_example(binary, workdir=None):
    """Run the TIP49 example that ships with SAINTexpress.

    The release includes the dataset the authors distribute with it, so this
    checks the program we are about to use against the output the published
    implementation produces on known data.
    """
    example = (Path(binary).resolve().parents[2]
               / "example input files" / "TIP49")
    if not example.is_dir():
        raise FileNotFoundError(f"no TIP49 example next to {binary}; expected {example}")

    temporary = workdir is None
    workdir = Path(workdir or tempfile.mkdtemp(prefix="coip_saint_tip49_"))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        for name in ("inter", "prey", "bait"):
            shutil.copy(example / f"{name}.dat", workdir / f"{name}.dat")
        done = subprocess.run([str(Path(binary).resolve()),
                               "inter.dat", "prey.dat", "bait.dat"],
                              cwd=workdir, capture_output=True, text=True)
        if done.returncode != 0:
            raise RuntimeError(f"the TIP49 example failed:\n{done.stderr}")
        return pd.read_csv(workdir / "list.txt", sep="\t")
    finally:
        if temporary:
            shutil.rmtree(workdir, ignore_errors=True)


# --- quality control on the counts ------------------------------------------

def count_summary(counts, condition_of):
    """Per sample: total spectra, proteins with at least one spectrum, depth.

    SAINTexpress compares a bait against its controls on the count scale, so a
    control that was simply sequenced less deeply looks like a cleaner control.
    These are the numbers needed to judge that.
    """
    return pd.DataFrame({
        "sample": counts.columns,
        "condition": [condition_of[s] for s in counts.columns],
        "total_spectra": counts.sum().astype(int).to_numpy(),
        "proteins_with_spectra": (counts > 0).sum().to_numpy(),
        "median_count_where_detected": [
            float(counts[s][counts[s] > 0].median()) for s in counts.columns],
        "max_count": counts.max().astype(int).to_numpy(),
    })


def validate_counts(counts):
    """Check the count matrix before it goes anywhere near SAINTexpress."""
    checks = []

    def check(name, passed, detail=""):
        checks.append({"check": name, "result": "pass" if passed else "FAIL",
                       "detail": detail})
        if not passed:
            raise ValueError(f"{name} , {detail}")

    values = counts.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    check("spectral counts are whole numbers",
          bool(np.all(finite == np.round(finite))))
    check("spectral counts are not negative", bool(np.all(finite >= 0)))
    check("no sample is empty", bool((counts.sum() > 0).all()))
    missing = int(np.isnan(values).sum())
    checks.append({"check": "missing versus zero counts", "result": "note",
                   "detail": f"{missing} missing cells, "
                             f"{int((values == 0).sum())} explicit zeros"})
    return pd.DataFrame(checks)


def tidy_scores(result, prey_mapping, baits):
    """Reshape SAINTexpress output to one row per protein group.

    SAINTexpress writes one row per bait-prey pair. This puts the baits side by
    side, indexed like the original workbook, so Analysis B can be joined to
    Analysis A without guessing which row is which.
    """
    key = prey_mapping.set_index("prey_id")["row"]
    wide = {}
    for bait in baits:
        panel = result[result["Bait"] == bait].copy()
        unknown = set(panel["Prey"]) - set(key.index)
        if unknown:
            raise ValueError(f"SAINTexpress returned unknown preys: {sorted(unknown)[:5]}")
        panel.index = key.reindex(panel["Prey"]).to_numpy()
        for source, name in (("AvgSpec", "avg_spec"), ("SpecSum", "spec_sum"),
                             ("Spec", "spec"), ("ctrlCounts", "control_counts"),
                             ("SaintScore", "saint_score"), ("AvgP", "avg_p"),
                             ("MaxP", "max_p"), ("FoldChange", "fold_change"),
                             ("BFDR", "bfdr")):
            wide[f"{name}_{bait}"] = panel[source]
    return pd.DataFrame(wide).reindex(prey_mapping["row"].to_numpy())
