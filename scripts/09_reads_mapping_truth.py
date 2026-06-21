#!/usr/bin/env python
"""Derive a species-level truth set from a CAMI reads_mapping.tsv (T10 / HMP).

The CAMI per_bodysite tarballs ship reads + a per-read truth map
(`reads/reads_mapping.tsv.gz`: anonymous_read_id, genome_id, tax_id, read_id)
but no species-level gold-standard .profile (unlike marine/strain). The source
tax_id is often strain-level, while Kraken2 reports species (rank S), so a naive
taxid set would under-count true positives. This rolls each distinct source
tax_id up to its species ancestor using the SAME nodes.dmp the classifier DB
uses, so truth and prediction share one taxonomy.

Usage:
    python scripts/09_reads_mapping_truth.py --mapping reads_mapping.tsv.gz \\
        --nodes databases/kraken2/k2_standard_08gb/nodes.dmp --out truth.txt
"""
import argparse
import gzip
from pathlib import Path


def load_nodes(path: Path) -> dict[str, tuple[str, str]]:
    """Parse nodes.dmp -> {taxid: (parent_taxid, rank)}."""
    nodes: dict[str, tuple[str, str]] = {}
    with open(path) as fh:
        for line in fh:
            cols = [c.strip() for c in line.split("|")]
            if len(cols) >= 3:
                nodes[cols[0]] = (cols[1], cols[2])
    return nodes


def to_rank(taxid: str, nodes: dict[str, tuple[str, str]], rank: str = "species") -> str | None:
    """Walk parent links to the ancestor at ``rank``; None if none/unknown.
    Stops at the root (parent == self) to avoid an infinite loop on a bad dump.

    Note: this only walks UP, so a taxid already coarser than ``rank`` (e.g. a
    genus taxid with rank='species') yields None — by design. CAMI 97%-OTU truth
    is often genus-resolved, so HMP body sites are scored with rank='genus'."""
    seen = set()
    cur = taxid
    while cur in nodes and cur not in seen:
        seen.add(cur)
        parent, cur_rank = nodes[cur]
        if cur_rank == rank:
            return cur
        if parent == cur:
            break
        cur = parent
    return None


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def reads_mapping_truth(mapping: Path, nodes: dict[str, tuple[str, str]], rank: str = "species") -> set[str]:
    raw: set[str] = set()
    with _open(mapping) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) >= 3:
                raw.add(cols[2].strip())
    return {r for t in raw if (r := to_rank(t, nodes, rank))}


def main() -> None:
    ap = argparse.ArgumentParser(description="CAMI reads_mapping.tsv -> truth set at a rank.")
    ap.add_argument("--mapping", required=True, type=Path)
    ap.add_argument("--nodes", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--rank", default="species", help="taxonomic rank for the truth set (species/genus/...)")
    args = ap.parse_args()
    nodes = load_nodes(args.nodes)
    truth = reads_mapping_truth(args.mapping, nodes, args.rank)
    args.out.write_text("\n".join(sorted(truth, key=int)) + "\n")
    print(f"{len(truth)} {args.rank} taxa -> {args.out}")


def _selfcheck() -> None:
    # tree: root1 > genus7 > species9 > {strain99, no-rank88}
    nodes = {"1": ("1", "no rank"), "7": ("1", "genus"), "9": ("7", "species"),
             "99": ("9", "strain"), "88": ("9", "no rank")}
    assert to_rank("99", nodes) == "9"            # strain -> species
    assert to_rank("88", nodes) == "9"            # dedups to same species
    assert to_rank("9", nodes) == "9"
    assert to_rank("7", nodes) is None            # genus has no species ancestor (walk up only)
    assert to_rank("404", nodes) is None          # unknown taxid
    assert to_rank("99", nodes, "genus") == "7"   # strain -> genus
    assert to_rank("9", nodes, "genus") == "7"    # species -> genus (the HMP case)
    print("selfcheck OK")


if __name__ == "__main__":
    import sys
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
