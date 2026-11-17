"""Fetch the two proteomes a species assignment has to be checked against.

The FR002 workbook names an organism for every accession in its FASTA_HEADERS
column, which is the search database's own annotation. That settles which
organism an accession belongs to, but not whether the peptides behind a
protein group could equally have come from the other species: human and
parasite share deeply conserved proteins, and a peptide that exists in both
cannot assign a protein group to either.

Answering that needs sequences, so this downloads:

  human     UniProtKB reviewed, taxon 9606 , the canonical human proteins
  parasite  UniProtKB taxon 5833 , P. falciparum including the 3D7 isolate
            and the strain-specific entries, reviewed and unreviewed
  retired   the newest archived version of any accession in the workbooks
            given on the command line that neither current set still carries

The retired set matters more than it sounds. A third of this experiment's
parasite accessions have been removed from UniProtKB since the search was
run, and without their sequences their peptides look as though they occur in
no parasite protein, which turns a conserved parasite protein into a false
human assignment.

Unreviewed human entries are left out; they are mostly fragments and isoform
variants of the reviewed sequences, and including the full 260k-entry set
would not change which peptides are shared between the two species.

    python build_species_proteomes.py [workbook.xlsx ...]

Prints the UniProt release and the SHA-256 of each file it writes. Re-running
it against a later release may change a few sequences, so record what the
analysis used.
"""

from __future__ import annotations

import gzip
import hashlib
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNISAVE = "https://rest.uniprot.org/unisave"
SETS = {
    "uniprot_human_reviewed.fasta.gz": "taxonomy_id:9606 AND reviewed:true",
    "uniprot_pfalciparum_5833.fasta.gz": "taxonomy_id:5833",
}


def get(url, attempts=5):
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=600) as response:
                return response.read(), response.headers
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503) or attempt == attempts - 1:
                raise
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"gave up on {url}")


def accessions_of(workbook):
    """Every accession named in a ProteoFlux Summary sheet."""
    import pandas as pd

    groups = pd.read_excel(workbook, sheet_name="Summary")["UNIPROT_AC"]
    return {a for group in groups for a in str(group).split(";")}


def entries_of(text):
    entries, header, parts = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                entries.append((header, "".join(parts)))
            header, parts = line, []
        elif header:
            parts.append(line.strip())
    if header:
        entries.append((header, "".join(parts)))
    return entries


def retired(wanted, present, directory,
            name="uniprot_retired_accessions.fasta.gz"):
    """Archived sequences for accessions no longer in the current releases."""
    missing = sorted(wanted - present)
    print(f"\n{len(missing)} accessions are not in the current releases")
    written, unresolved = 0, []
    with gzip.open(directory / name, "wt") as out:
        for n, accession in enumerate(missing, 1):
            # UniSave returns every archived version, newest first.
            payload, _ = get(f"{UNISAVE}/{accession}?format=fasta")
            entries = entries_of(payload.decode())
            if entries:
                header, sequence = entries[0]
                out.write(f">sp|{accession}|archived {header[1:]}\n")
                for i in range(0, len(sequence), 60):
                    out.write(sequence[i:i + 60] + "\n")
                written += 1
            else:
                unresolved.append(accession)
            if n % 25 == 0:
                print(f"  UniSave {n}/{len(missing)}", flush=True)
            time.sleep(0.1)
    print(f"{name}: {written} archived sequences")
    if unresolved:
        print("  no sequence at all for:", ", ".join(unresolved))
    print("  SHA-256:",
          hashlib.sha256((directory / name).read_bytes()).hexdigest())


def main(directory=".", workbooks=()):
    directory = Path(directory)
    present = set()
    for name, query in SETS.items():
        parameters = urllib.parse.urlencode(
            {"query": query, "format": "fasta", "compressed": "true",
             "includeIsoform": "false"})
        payload, headers = get(f"{STREAM}?{parameters}")
        path = directory / name
        path.write_bytes(payload)
        entries = gzip.decompress(payload).decode().count(">")
        print(f"{name}: {entries:,} entries, "
              f"UniProt release {headers.get('X-UniProt-Release')} "
              f"({headers.get('X-UniProt-Release-Date')})")
        print("  SHA-256:", hashlib.sha256(payload).hexdigest())
        present |= {line.split("|")[1] for line
                    in gzip.decompress(payload).decode().splitlines()
                    if line.startswith(">")}

    wanted = set()
    for workbook in workbooks:
        wanted |= accessions_of(workbook)
    if wanted:
        retired(wanted, present, directory)


if __name__ == "__main__":
    main(Path(__file__).parent, sys.argv[1:])
