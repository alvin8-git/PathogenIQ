#!/usr/bin/env python3
"""
Download Tier-1 clinical pathogen genomes from NCBI RefSeq and build
a sourmash SBT database for PathogenIQ.

Uses concurrent.futures for parallel downloads + ncbi-genome-download
for the actual genome fetching.

Runtime: ~30-90 min depending on connection and number of genomes
Output:  databases/tier1/tier1_pathogens.sbt.zip

Prerequisites:
    pip install ncbi-genome-download sourmash requests

Usage:
    python scripts/02_download_tier1_db.py [--threads 16] [--outdir databases/tier1]
"""

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


# ---------------------------------------------------------------------------
# Tier-1 pathogen list — organised by clinical priority
# Drawn from: CZID pathogen list, WHO Priority Pathogens 2024,
#             CDC Select Agents, common clinical metagenomics panels
#
# Format: "NCBI taxonomy name": "RefSeq representative accession"
# ---------------------------------------------------------------------------
TIER1_PATHOGENS = {
    # ── Gram-negative ESKAPE / WHO Critical Priority ──────────────────────
    "Acinetobacter baumannii":          "GCF_000015425.1",
    "Klebsiella pneumoniae":            "GCF_000240185.1",
    "Pseudomonas aeruginosa":           "GCF_000006765.1",
    "Escherichia coli":                 "GCF_000005845.2",
    "Enterobacter cloacae":             "GCF_000025565.1",
    "Salmonella enterica":              "GCF_000006945.2",
    "Haemophilus influenzae":           "GCF_000001405.40",  # note: this is human, use below
    # fix — correct H. influenzae
    "Haemophilus influenzae Rd KW20":   "GCF_000001469.1",
    "Neisseria meningitidis":           "GCF_000008805.1",
    "Neisseria gonorrhoeae":            "GCF_000006845.1",
    "Burkholderia pseudomallei":        "GCF_000011545.1",  # CDC Select Agent
    "Yersinia pestis":                  "GCF_000009065.1",  # CDC Select Agent
    "Francisella tularensis":           "GCF_000008985.1",  # CDC Select Agent

    # ── Gram-positive ─────────────────────────────────────────────────────
    "Staphylococcus aureus":            "GCF_000013425.1",
    "Enterococcus faecalis":            "GCF_000007785.1",
    "Enterococcus faecium":             "GCF_000174395.2",
    "Streptococcus pneumoniae":         "GCF_000007045.1",
    "Streptococcus pyogenes":           "GCF_000006785.2",
    "Listeria monocytogenes":           "GCF_000196035.1",
    "Bacillus anthracis":               "GCF_000008445.1",  # CDC Select Agent
    "Clostridioides difficile":         "GCF_000009205.2",
    "Mycobacterium tuberculosis":       "GCF_000195955.2",
    "Mycobacterium abscessus":          "GCF_000069185.1",

    # ── Anaerobes ─────────────────────────────────────────────────────────
    "Bacteroides fragilis":             "GCF_000025985.1",
    "Fusobacterium nucleatum":          "GCF_000007325.1",

    # ── Spirochetes ───────────────────────────────────────────────────────
    "Borrelia burgdorferi":             "GCF_000008685.2",
    "Treponema pallidum":               "GCF_000008605.1",
    "Leptospira interrogans":           "GCF_000007685.1",

    # ── Atypicals ─────────────────────────────────────────────────────────
    "Legionella pneumophila":           "GCF_000048645.1",
    "Chlamydia trachomatis":            "GCF_000008585.1",
    "Mycoplasma pneumoniae":            "GCF_000027345.1",

    # ── DNA Viruses ───────────────────────────────────────────────────────
    "Human herpesvirus 1":              "GCF_000859985.2",   # HSV-1
    "Human herpesvirus 2":              "GCF_000861065.1",   # HSV-2
    "Human herpesvirus 3":              "GCF_000858285.1",   # VZV
    "Human herpesvirus 4":              "GCF_000862125.1",   # EBV
    "Human herpesvirus 5":              "GCF_000845245.1",   # CMV
    "Human herpesvirus 6A":             "GCF_000848125.1",
    "Human herpesvirus 8":              "GCF_000863945.1",   # KSHV
    "Adenovirus C":                     "GCF_000858705.1",
    "BK polyomavirus":                  "GCF_000864765.1",
    "JC polyomavirus":                  "GCF_000864105.1",
    "Hepatitis B virus":                "GCF_000861765.1",   # NC_003977 ayw strain
    "Human papillomavirus 16":          "GCF_000863645.1",

    # ── RNA Viruses ───────────────────────────────────────────────────────
    "SARS-CoV-2":                       "GCF_009858895.2",
    "Influenza A virus H1N1":           "GCF_001343785.1",
    # TODO: add correct GCF accessions for these RNA viruses once verified:
    #   "Influenza B virus"  — was GCF_000864765.1 (duplicate of BK polyomavirus)
    #   "Hepatitis C virus"  — was GCF_000864765.1 (duplicate of BK polyomavirus)
    #   "HIV-1"              — was GCF_000864765.1 (duplicate of BK polyomavirus)
    #   "Dengue virus"       — was GCF_000862125.1 (duplicate of EBV)
    #   "West Nile virus"    — was GCF_000864105.1 (duplicate of JC polyomavirus)
    "Rabies lyssavirus":                "GCF_000859085.1",
    "Enterovirus A":                    "GCF_000859065.1",   # EV-A71

    # ── Fungi ─────────────────────────────────────────────────────────────
    "Candida albicans":                 "GCF_000182965.3",
    "Candida auris":                    "GCF_002775015.1",
    "Aspergillus fumigatus":            "GCF_000002655.1",
    "Cryptococcus neoformans":          "GCF_000149245.2",
    "Pneumocystis jirovecii":           "GCF_001477535.1",
    "Histoplasma capsulatum":           "GCF_000006445.1",
    "Coccidioides immitis":             "GCF_000149335.2",
    "Mucor circinelloides":             "GCF_000401635.1",
    "Saccharomyces cerevisiae":         "GCF_000146045.2",   # ZymoBIOMICS

    # ── Parasites ─────────────────────────────────────────────────────────
    "Plasmodium falciparum":            "GCF_000002765.5",
    "Toxoplasma gondii":                "GCF_000006565.2",
    "Cryptosporidium parvum":           "GCF_000165345.1",
    "Giardia lamblia":                  "GCF_000002435.1",
    "Trypanosoma brucei":               "GCF_000002445.2",
    "Leishmania major":                 "GCF_000002725.1",
}

# Accessions with known issues — skip these
_SKIP = {
    "GCF_000001405.40",   # GRCh38 (human genome), not H. influenzae
    # Not available via ncbi-genome-download (refseq + genbank both return 0 results).
    # Replace with updated accessions when found:
    "GCF_000001469.1",    # H. influenzae Rd KW20 — try GCF_000001469.2
    "GCF_000863945.1",    # HHV-8/KSHV
    "GCF_000149245.2",    # Cryptococcus neoformans JEC21
    "GCF_000006445.1",    # Histoplasma capsulatum
    "GCF_000002725.1",    # Leishmania major Friedlin
    "GCF_000401635.1",    # Mucor circinelloides CBS277.49
}

VALID_ACCESSIONS = {
    name: acc for name, acc in TIER1_PATHOGENS.items()
    if acc not in _SKIP
}

def download_genome(name: str, accession: str, genome_dir: Path, retries: int = 3) -> tuple[str, bool, str]:
    """Download a single genome; tries refseq then genbank."""
    if list(genome_dir.rglob(f"{accession}*.fna.gz")):
        return name, True, f"already present ({accession})"
    for section in ("refseq", "genbank"):
        cmd = [
            "ncbi-genome-download",
            "-A", accession,
            "-F", "fasta",
            "-o", str(genome_dir),
            "-r", str(retries),
            "-s", section,
            "-l", "all",
            "all",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if list(genome_dir.rglob(f"{accession}*.fna.gz")):
            return name, True, f"downloaded from {section} ({accession})"
    return name, False, f"FAILED ({accession}): {result.stderr.strip()[:120]}"

def sketch_genome(fna_path: Path, sig_dir: Path) -> tuple[Path, bool]:
    """Sketch a single genome with sourmash. Returns (path, success)."""
    sig_out = sig_dir / (fna_path.stem + ".sig")
    if sig_out.exists():
        return fna_path, True
    cmd = [
        "sourmash", "sketch", "dna",
        "-p", "k=31,scaled=1000",
        str(fna_path),
        "-o", str(sig_out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return fna_path, result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="Download Tier-1 pathogen DB for PathogenIQ")
    parser.add_argument("--threads", type=int, default=8, help="Parallel download threads")
    parser.add_argument("--outdir", default="databases/tier1", help="Output directory")
    parser.add_argument("--zymo-only", action="store_true", help="Download ZymoBIOMICS 10 organisms only")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    genome_dir = outdir / "genomes"
    sig_dir = outdir / "sigs"
    db_out = outdir / "tier1_pathogens.sbt.zip"
    genome_dir.mkdir(parents=True, exist_ok=True)
    sig_dir.mkdir(parents=True, exist_ok=True)

    targets = VALID_ACCESSIONS
    if args.zymo_only:
        zymo_names = {
            "Pseudomonas aeruginosa", "Escherichia coli", "Salmonella enterica",
            "Lactobacillus fermentum ATCC14931", "Enterococcus faecalis",
            "Staphylococcus aureus", "Listeria monocytogenes",
            "Bacillus subtilis 168", "Saccharomyces cerevisiae", "Cryptococcus neoformans",
        }
        targets = {k: v for k, v in VALID_ACCESSIONS.items()
                   if any(n in k for n in zymo_names)}

    print(f"Downloading {len(targets)} genomes with {min(args.threads, 4)} threads...")
    print(f"Output directory: {outdir}\n")

    # ── Step 1: Per-accession downloads (isolated failures, refseq→genbank fallback) ──
    failed = []
    dl_workers = min(args.threads, 4)  # >4 triggers NCBI rate-limiting
    with ThreadPoolExecutor(max_workers=dl_workers) as executor:
        futures = {
            executor.submit(download_genome, name, acc, genome_dir): name
            for name, acc in targets.items()
        }
        for i, future in enumerate(as_completed(futures), 1):
            name, success, msg = future.result()
            status = "✓" if success else "✗"
            print(f"  [{i:>3}/{len(targets)}] {status} {name}: {msg}")
            if not success:
                failed.append(name)

    if failed:
        print(f"\nWARNING: {len(failed)} downloads failed:")
        for f in failed:
            print(f"  - {f}")

    # ── Step 2: Decompress genomes ─────────────────────────────────────────
    print("\nDecompressing genomes...")
    gz_files = list(genome_dir.rglob("*.fna.gz"))
    for gz in gz_files:
        fna = gz.parent / gz.stem  # removes .gz, keeps full path
        if not fna.exists():
            result = subprocess.run(["gunzip", "-k", str(gz)], capture_output=True)
            if result.returncode != 0:
                print(f"  WARNING: skipping corrupt/non-gzip file: {gz.name}")

    # ── Step 3: Parallel sketching ─────────────────────────────────────────
    fna_files = list(genome_dir.rglob("*.fna"))
    print(f"\nSketching {len(fna_files)} genome files...")

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(sketch_genome, f, sig_dir): f for f in fna_files}
        for i, future in enumerate(as_completed(futures), 1):
            fna_path, success = future.result()
            status = "✓" if success else "✗"
            print(f"  [{i:>3}/{len(fna_files)}] {status} {fna_path.name}")

    # ── Step 4: Build SBT ──────────────────────────────────────────────────
    sig_files = list(sig_dir.glob("*.sig"))
    print(f"\nBuilding sourmash SBT from {len(sig_files)} signatures...")
    subprocess.run(
        ["sourmash", "index", "--ksize", "31", str(db_out)] + [str(s) for s in sig_files],
        check=True,
    )

    print(f"\n{'='*60}")
    print(f"Tier-1 DB ready: {db_out}")
    print(f"Genomes included: {len(sig_files)}")
    print(f"\nSet environment variable:")
    print(f"  export TIER1_DB=$(pwd)/{db_out}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
