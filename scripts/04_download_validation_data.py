#!/usr/bin/env python
"""Download validation / NTC-background data for the NTC-aware grading work (T8).

Datasets
--------
- salter-ntc : Salter et al. 2014 metagenomic reagent/kitome controls
               (ENA project ERP006808). The ultrapure-water negative controls
               are the source for the Tier-2 pooled background_default.tsv
               (feed them to scripts/build_background_default.py afterwards).
- salter-16s : Salter et al. 2014 16S controls (ENA ERP006737) — amplicon,
               usually not needed for shotgun grading; included for completeness.
- cami       : CAMI II gold-standard multi-community truth sets. The full set is
               large and spread across many records; this script does not guess
               record URLs. Pass the specific file URLs you want with --url, using
               the portals documented below.

CAMI entry points (pick specific records, then pass their file URLs via --url):
  - https://cami-challenge.org/                (schedule / data pages)
  - https://zenodo.org/communities/cami        (results + gold standards)
  - GigaDB DOI 10.5524/100344                  (benchmark datasets)

ENA accessions are resolved to FASTQ URLs automatically via the ENA filereport
API, with md5 verification and resumable downloads.

Examples
--------
  # list what the Salter NTC project contains without downloading
  python scripts/04_download_validation_data.py --dataset salter-ntc --dry-run

  # download the Salter reagent/kitome controls
  python scripts/04_download_validation_data.py --dataset salter-ntc --out databases/validation

  # download specific CAMI files you picked from the portals above
  python scripts/04_download_validation_data.py --dataset cami \\
      --url https://zenodo.org/record/<ID>/files/<file>.tar.gz

Stdlib only — no extra dependencies. Resumes partial downloads (.part files).
"""
import argparse
import hashlib
import urllib.parse
import urllib.request
from pathlib import Path

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

DATASETS = {
    "salter-ntc": {
        "kind": "ena", "accession": "ERP006808",
        "desc": "Salter 2014 metagenomic reagent/kitome controls (NTC background source)",
    },
    "salter-16s": {
        "kind": "ena", "accession": "ERP006737",
        "desc": "Salter 2014 16S controls (amplicon)",
    },
    "cami": {
        "kind": "urls", "urls": [],
        "desc": "CAMI II truth sets — pass file URLs via --url (see portals in the module docstring)",
    },
    # Non-spiked negative-control shotgun datasets for a cleaner, broader Tier-2
    # kitome background (pool with scripts/11_pool_blank_background.py). Only the
    # blank runs are downloaded — title_filter selects them from projects that
    # also contain real samples.
    "blanks-hunt": {
        "kind": "ena", "accession": "PRJEB66439", "title_filter": "Blank_sample_DNA",
        "desc": "HUNT One Health reagent/extraction blanks (47 runs, RANDOM shotgun, ~1.2M reads each)",
    },
    "blanks-qiita": {
        "kind": "ena", "accession": "PRJEB56784", "title_filter": "BLANK",
        "desc": "Qiita 14332 nucleic-acid pipeline reagent blanks (EtOH/IPA/swab/buffer, PCR-WGS, shallow)",
    },
    # Air / bioaerosol surveillance prep (Plan 5). PRJNA1228129 = aircraft-cabin
    # air filters + masks; the only shortlisted air dataset with BOTH negative
    # controls and positive spike-ins. See docs/air-datasets-2026-06-22.md.
    "air-aircraft": {
        "kind": "ena", "accession": "PRJNA1228129",
        "desc": "Aircraft-cabin air metagenomes — all 55 runs (air filters/masks + controls + spike-ins)",
    },
    "air-aircraft-ntc": {
        "kind": "ena", "accession": "PRJNA1228129", "title_filter": "Control",
        "desc": "Aircraft-cabin air NEGATIVE CONTROLS (air kitome NTC for build_background)",
    },
    "air-aircraft-spike": {
        "kind": "ena", "accession": "PRJNA1228129", "title_filter": "Spike",
        "desc": "Aircraft-cabin air SPIKE-INS (known positives for air detection precision/recall)",
    },
}


def _filename_from_url(url: str) -> str:
    """Derive a sensible filename from a download URL. Handles the Zenodo pattern
    .../files/<filename>/content, where the real name is the second-to-last
    path segment rather than the literal 'content'."""
    parts = urllib.parse.urlparse(url).path.rstrip("/").split("/")
    name = parts[-1] if parts else ""
    if name in ("content", "download") and len(parts) >= 2:
        name = parts[-2]
    return name or "download"


def _parse_filereport(text: str, title_filter: str | None = None) -> list[tuple[str, str, str]]:
    """Parse an ENA filereport TSV into (run_accession, https_url, md5) rows.

    Columns: run_accession, fastq_ftp, fastq_md5, sample_title. fastq_ftp /
    fastq_md5 are ';'-separated and paired (paired-end runs have two files). ENA
    serves ftp.sra.ebi.ac.uk/... paths over HTTPS too, so we prefix https:// when
    no scheme is present. ``title_filter`` keeps only rows whose sample_title
    contains that substring (case-insensitive) — used to pull just the blank runs
    from a project that also holds real samples.
    """
    files: list[tuple[str, str, str]] = []
    lines = [ln for ln in text.splitlines() if ln.strip()]
    for ln in lines[1:]:   # skip header
        cols = ln.split("\t")
        if len(cols) < 3:
            continue
        run, ftp, md5 = cols[0], cols[1], cols[2]
        title = cols[3] if len(cols) > 3 else ""
        if title_filter and title_filter.lower() not in title.lower():
            continue
        urls, md5s = ftp.split(";"), md5.split(";")
        for u, m in zip(urls, md5s):
            if not u:
                continue
            full = u if u.startswith(("http://", "https://", "ftp://")) else "https://" + u
            files.append((run, full, m))
    return files


def ena_runs(accession: str, title_filter: str | None = None) -> list[tuple[str, str, str]]:
    params = urllib.parse.urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_md5,sample_title",
        "format": "tsv",
    })
    with urllib.request.urlopen(f"{ENA_FILEREPORT}?{params}", timeout=60) as r:
        return _parse_filereport(r.read().decode(), title_filter)


def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def download(url: str, dest: Path, expected_md5: str | None = None) -> None:
    dest = Path(dest)
    if dest.exists() and expected_md5 and md5sum(dest) == expected_md5:
        print(f"  = {dest.name} (md5 ok, skip)")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    pos = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(url)
    if pos:
        req.add_header("Range", f"bytes={pos}-")
    with urllib.request.urlopen(req, timeout=300) as r:
        # If the server ignored our Range (status 200), restart from scratch.
        resuming = pos and getattr(r, "status", 200) == 206
        mode = "ab" if resuming else "wb"
        print(f"  > {dest.name}" + (f" (resume @ {pos})" if resuming else ""))
        with open(part, mode) as f:
            while True:
                blk = r.read(1 << 20)
                if not blk:
                    break
                f.write(blk)
    if expected_md5:
        got = md5sum(part)
        if got != expected_md5:
            raise SystemExit(f"md5 mismatch for {dest.name}: {got} != {expected_md5}")
    part.rename(dest)
    print(f"  + {dest.name}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Download validation/NTC data (T8).")
    ap.add_argument("--dataset", choices=list(DATASETS) + ["all"], default="salter-ntc")
    ap.add_argument("--out", default="databases/validation", help="output directory")
    ap.add_argument("--url", action="append", default=[],
                    help="extra direct URL(s) to download (e.g. CAMI record files)")
    ap.add_argument("--dry-run", action="store_true", help="list files, download nothing")
    ap.add_argument("--limit", type=int, default=0,
                    help="for ENA datasets, download only the first N runs (0 = all)")
    args = ap.parse_args()

    names = list(DATASETS) if args.dataset == "all" else [args.dataset]
    for name in names:
        d = DATASETS[name]
        outdir = Path(args.out) / name
        print(f"[{name}] {d['desc']}")
        if d["kind"] == "ena":
            files = ena_runs(d["accession"], d.get("title_filter"))
            if args.limit:
                keep = list(dict.fromkeys(run for run, _u, _m in files))[:args.limit]
                files = [(r, u, m) for r, u, m in files if r in set(keep)]
            tf = f" matching '{d['title_filter']}'" if d.get("title_filter") else ""
            lim = f" (limit {args.limit} runs)" if args.limit else ""
            print(f"  {len(files)} file(s) in {d['accession']}{tf}{lim}")
            for _run, url, md5 in files:
                fn = _filename_from_url(url)
                if args.dry_run:
                    print(f"  would download {url}")
                    continue
                download(url, outdir / fn, md5)
        else:
            urls = list(d.get("urls", [])) + args.url
            if not urls:
                print("  no URLs — pass CAMI file URLs with --url (see module docstring)")
            for url in urls:
                fn = _filename_from_url(url)
                if args.dry_run:
                    print(f"  would download {url}")
                    continue
                download(url, outdir / fn)

    # extra --url for an ENA dataset selection (downloaded alongside)
    if args.url and DATASETS.get(args.dataset, {}).get("kind") != "urls":
        for url in args.url:
            fn = _filename_from_url(url)
            if args.dry_run:
                print(f"would download {url}")
                continue
            download(url, Path(args.out) / "extra" / fn)


if __name__ == "__main__":
    main()
