# PathogenIQ TODO — Next Session

## Critical Bug: Pipeline Produces Broken Output

**Status:** Root cause identified, fix not yet implemented.

**Problem:** Running the full pipeline on Zymo data produces reports with:
- Empty organism names
- NaN abundances and CI bounds
- INT64_MIN read counts (`-9223372036854775808`)
- All Grade X

**Root Cause (in `pathogeniq/sketch.py`):**
When `sourmash search` queries an SBT (Sequence Bloom Tree) index, the CSV output has:
- `name` column: **empty** (signatures were built without `--name`)
- `filename` column: the SBT zip file path, not individual genome paths
- `md5` column: valid md5 hash of the matched signature

This causes `SketchHit.name=""` and `SketchHit.genome_path=SBT.zip`, so:
1. `align.py` passes empty names to the report
2. `align.py` runs `minimap2` against a `.zip` file → empty PAFs
3. EM receives zero-alignment matrix → NaN/underflow values

**Fix Required:**
Build an md5-to-genome mapping from the signature files in `databases/tier1/sigs/`.
Each `.sig` file's md5 (from `sourmash sig describe`) maps to:
- Genome name = sig filename stem (e.g. `GCF_000005845.2_ASM584v2_genomic`)
- Genome path = `databases/tier1/genomes/{stem}.fna`

In `sketch.py`, after parsing CSV rows, look up each row's `md5` in the mapping
to get the correct `name` and `genome_path`.

## Already Completed This Session

- [x] README.md updated with Plan 2 features (AMR, PDF, contaminants, pipeline diagram)
- [x] README.md committed and pushed to GitHub
- [x] Version.md and Documentation.md created in previous session
- [x] Zymo pipeline output inspected — bug identified

## Next Steps (Priority Order)

1. **Fix `pathogeniq/sketch.py`** — md5-to-genome mapping for SBT search results
2. **Re-run Zymo pipeline** to verify fix
3. **Run unit tests** (`pytest tests/ -v`) to ensure no regressions
4. **Commit and push** the fix
