#!/usr/bin/env python
"""Build a custom Kraken2 DB from the tier1 genomes (for the T10 benchmark).

Each genome's species is read from name_map.json (keyed by the GCF accession in
the filename), resolved to an NCBI taxid via the downloaded taxonomy names.dmp,
and every sequence header is tagged with |kraken:taxid|TAXID so Kraken2 can build
without the huge accession2taxid maps. Then the library is added and the DB built.

Prereq:
  conda run -n pathogeniq kraken2-build --download-taxonomy --skip-maps \\
      --db databases/kraken2_tier1

Usage (in the pathogeniq env, so kraken2-build is on PATH):
  conda run -n pathogeniq python scripts/07_build_kraken_db.py \\
      --db databases/kraken2_tier1 --genomes databases/tier1/genomes \\
      --name-map databases/tier1/name_map.json --threads 16
"""
import argparse
import json
import re
import subprocess
from pathlib import Path

_GCF = re.compile(r"(GC[AF]_\d+\.\d+)")


def load_name2taxid(names_dmp: Path) -> dict[str, int]:
    """Scientific-name (lowercased) -> taxid from names.dmp."""
    out: dict[str, int] = {}
    for line in names_dmp.read_text().splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) >= 4 and parts[3] == "scientific name":
            out[parts[1].lower()] = int(parts[0])
    return out


def resolve_taxid(species: str | None, name2taxid: dict[str, int]) -> int | None:
    if not species:
        return None
    t = name2taxid.get(species.lower())
    if t is not None:
        return t
    # fall back to the genus (first token) for strain-named species
    return name2taxid.get(species.split()[0].lower())


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a custom Kraken2 DB from tier1 genomes.")
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--genomes", required=True, type=Path)
    ap.add_argument("--name-map", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    name2taxid = load_name2taxid(args.db / "taxonomy" / "names.dmp")
    name_map = json.loads(args.name_map.read_text())   # GCF accession -> species name

    staging = args.db / "staging"
    staging.mkdir(parents=True, exist_ok=True)

    added = 0
    unmapped: list[tuple[str, str | None]] = []
    for fna in sorted(args.genomes.glob("*.fna")):
        m = _GCF.search(fna.name)
        if not m:
            continue
        species = name_map.get(m.group(1))
        taxid = resolve_taxid(species, name2taxid)
        if taxid is None:
            unmapped.append((fna.name, species))
            continue
        tagged = staging / fna.name
        with open(fna) as fin, open(tagged, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    seqid = line[1:].split()[0]
                    fout.write(f">{seqid}|kraken:taxid|{taxid}\n")
                else:
                    fout.write(line)
        subprocess.run(
            ["kraken2-build", "--add-to-library", str(tagged), "--db", str(args.db)],
            check=True, capture_output=True,
        )
        added += 1

    print(f"added {added} genomes to the library; {len(unmapped)} unmapped")
    for name, species in unmapped:
        print(f"  unmapped: {name} -> {species}")

    print("building DB (this is the slow step)...")
    subprocess.run(
        ["kraken2-build", "--build", "--db", str(args.db), "--threads", str(args.threads)],
        check=True,
    )
    print(f"DB built at {args.db}")


if __name__ == "__main__":
    main()
