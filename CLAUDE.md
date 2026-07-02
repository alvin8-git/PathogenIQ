# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install (editable, with dev deps)
pip install -e ".[dev]"

# Unit tests (fast, no external tools needed)
pytest tests/ -v

# Single test
pytest tests/test_report.py::test_grade_a -v

# Integration tests (requires databases + external tools)
pytest tests/integration/ -v -m integration

# Lint
ruff check pathogeniq/ tests/

# Type check
mypy pathogeniq/
```

## Architecture

PathogenIQ is a clinical metagenomics pipeline. Each stage is a separate module that wraps an external bioinformatics tool and returns a typed dataclass.

**Entrypoint:** `cli.py` is a `click` group; the `run` command orchestrates the pipeline. Required options: `--input` (FASTQ), `--output` (dir), `--db` (Tier-1 DB), `--host-ref`, `--specimen`. Notable extras: `--read-type short|long`, `--amr-db` (default `card`), `--ntc`/`--background`/`--no-background` (NTC tier), `--spike-taxon`/`--spike-copies`/`--sample-volume` (absolute quant), `--novelty`, `--viral`, `--assemble` (open-world arms), `--no-pdf`, `--sketch-threshold`, `--n-bootstrap`.

**Pipeline flow:**
```
QC (fastp/Chopper) → Host removal (BWA-MEM2/minimap2) → PhiX removal (minimap2)
  → [--novelty: Kraken2 broad DB] → Sketch screen (sourmash)
  → Targeted alignment (minimap2, +breadth-of-coverage) → EM abundance
  → AMR + virulence (ABRicate CARD/VFDB) → NTC background subtraction (background.py)
  → [--spike-*: absolute quant] → Grading (A/B/C/X)
  → [--assemble: MEGAHIT→MetaBAT2→CheckM→GTDB-Tk + pathogenicity triage]
  → [--viral: MEGAHIT→geNomad→CheckV]
  → Report (JSON + TSV + PDF + HTML)
```
The open-world arms (`--assemble`/`--viral`) run via `_discovery_arms()` in `cli.py` even when the targeted screen finds zero hits (a viral/novel-only sample still produces a report).

**Key data contracts between stages:**
- `PipelineConfig` (`config.py`) — config dataclass passed to every stage; `ReadType` (short/long), `SpecimenType` (blood/csf/bal/tissue/**air**), paths, tuning params
- `AlignResult` (`align.py`) — `alignment_matrix` + `organism_names` + `taxon_ids` + **`coverage`** (per-organism `CoverageStats`); bridges alignment → EM and feeds the breadth gate
- `EMResult` (`em.py`) → consumed by `build_entries()`/`write_report()`
- `BackgroundModel` (`background.py`) — per-taxon RPM rates + tier; `is_background()` (NB upper tail) zeroes out kitome, EXCEPT `is_dual_use()` taxa which are flagged not subtracted
- `ReportEntry` (`report.py`) — one row per targeted organism; lazy `grade` property reads `breadth_ratio`, `contaminant_risk`, tier, crossmap
- Open-world results: `MAG` (`assembly.py`), `ViralContig` (`viral.py`), `PathogenicityAssessment` (`pathogenicity.py`), `NoveltyResult` (`novelty.py`) — each graded by `grade_open_world()` (capped at B) and emitted as its own JSON block
- `AMRHit`/`VirulenceHit` (`amr.py`) — contig-based; `map_contigs_to_organisms()` aligns contigs back to the screened genomes so each hit attributes to a finding via `organism_matches` (co-attributed to every sibling within `_COMAP_MARGIN`=0.70 of best aligned-bp, so E. coli/Shigella shared genes show under each). `SpikeInfo` (`quantify.py`)
- Renderers: `report.py` (JSON+TSV), `pdf_report.py`, `html_report.py` share the same `ReportEntry` rows

**Evidence grading** (`report.py` + `contaminants.py` + `coverage.py` + `crossmap.py`):
- Grade A/B/C/X is a lazy property on `ReportEntry` (pure `grade(GradingInput)`)
- Gates: specimen `_MIN_READS` floor · CI ≤ 0.15 + **Tier-1 NTC** for A · `breadth_ratio` < 0.25 → X · cross-mapping artifact → X · contaminant demotion
- `contaminants.py::CONTAMINANT_PRIORS` is per-`SpecimenType` (AIR prior is deliberately NOT genus-broad on Pseudomonas/Acinetobacter so real pathogens aren't demoted)
- Reference-free hits (MAGs/viral) use `grade_open_world()` — genome completeness + contamination + marker signal, capped at B

**Read-type branching:**
- Short reads: fastp for QC, BWA-MEM2 for host removal
- Long reads: Chopper for QC, minimap2 for host removal
- Both paths produce identical intermediate files; downstream stages are read-type agnostic

**External tool dependencies:** core (fastp, Chopper, BWA-MEM2, minimap2, samtools, sourmash) via `environment.yml`. Optional/open-world: kraken2 + broad DB (novelty), abricate (AMR/VFDB/R4 markers), megahit/metabat2/checkm/gtdbtk (MAG arm), genomad/checkv (viral arm). geNomad/CheckV/abricate need `numpy<2` and live in **isolated conda envs** symlinked/wrapped onto `PATH` (`scripts/13_setup_viral_env.sh`) — never install them into the core env (it breaks SciPy). The MAG toolchain (metabat2/checkm/gtdbtk) uses the same discipline via `scripts/15_setup_mag_env.sh` (isolated `mag-bin` + `gtdbtk` envs, DBs on a data volume); `scripts/16_run_air_assemble.sh` runs the `--assemble` arm on the aircraft filters. Every external stage is **non-blocking**: a missing tool/DB skips that stage rather than failing the run.

## Test structure

- `tests/test_*.py` — unit tests; mock all subprocess/tool calls; run without any databases or external tools
- `tests/integration/test_zymo_validation.py` — end-to-end tests against ZymoBIOMICS D6300 standard; marked `@pytest.mark.integration`; require `HUMAN_REF` and `TIER1_DB` env vars pointing to built databases

## Database setup

Three databases must be built before running the pipeline or integration tests:

```bash
bash scripts/01_download_zymo_db.sh 16        # ZymoBIOMICS validation DB (~5 min)
python scripts/02_download_tier1_db.py --threads 16  # Tier-1 clinical DB, ~110 pathogen genomes (~30-90 min)
bash scripts/03_download_human_ref.sh 16       # GRCh38 + index (~90 min, needs ~64 GB RAM)
```

Output lands in `databases/`. See `scripts/README.md` for accession tables.

# context-mode — MANDATORY routing rules

You have context-mode MCP tools available. These rules are NOT optional — they protect your context window from flooding. A single unrouted command can dump 56 KB into context and waste the entire session.

## BLOCKED commands — do NOT attempt these

### curl / wget — BLOCKED
Any Bash command containing `curl` or `wget` is intercepted and replaced with an error message. Do NOT retry.
Instead use:
- `ctx_fetch_and_index(url, source)` to fetch and index web pages
- `ctx_execute(language: "javascript", code: "const r = await fetch(...)")` to run HTTP calls in sandbox

### Inline HTTP — BLOCKED
Any Bash command containing `fetch('http`, `requests.get(`, `requests.post(`, `http.get(`, or `http.request(` is intercepted and replaced with an error message. Do NOT retry with Bash.
Instead use:
- `ctx_execute(language, code)` to run HTTP calls in sandbox — only stdout enters context

### WebFetch — BLOCKED
WebFetch calls are denied entirely. The URL is extracted and you are told to use `ctx_fetch_and_index` instead.
Instead use:
- `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` to query the indexed content

## REDIRECTED tools — use sandbox equivalents

### Bash (>20 lines output)
Bash is ONLY for: `git`, `mkdir`, `rm`, `mv`, `cd`, `ls`, `npm install`, `pip install`, and other short-output commands.
For everything else, use:
- `ctx_batch_execute(commands, queries)` — run multiple commands + search in ONE call
- `ctx_execute(language: "shell", code: "...")` — run in sandbox, only stdout enters context

### Read (for analysis)
If you are reading a file to **Edit** it → Read is correct (Edit needs content in context).
If you are reading to **analyze, explore, or summarize** → use `ctx_execute_file(path, language, code)` instead. Only your printed summary enters context. The raw file content stays in the sandbox.

### Grep (large results)
Grep results can flood context. Use `ctx_execute(language: "shell", code: "grep ...")` to run searches in sandbox. Only your printed summary enters context.

## Tool selection hierarchy

1. **GATHER**: `ctx_batch_execute(commands, queries)` — Primary tool. Runs all commands, auto-indexes output, returns search results. ONE call replaces 30+ individual calls.
2. **FOLLOW-UP**: `ctx_search(queries: ["q1", "q2", ...])` — Query indexed content. Pass ALL questions as array in ONE call.
3. **PROCESSING**: `ctx_execute(language, code)` | `ctx_execute_file(path, language, code)` — Sandbox execution. Only stdout enters context.
4. **WEB**: `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` — Fetch, chunk, index, query. Raw HTML never enters context.
5. **INDEX**: `ctx_index(content, source)` — Store content in FTS5 knowledge base for later search.

## Subagent routing

When spawning subagents (Agent/Task tool), the routing block is automatically injected into their prompt. Bash-type subagents are upgraded to general-purpose so they have access to MCP tools. You do NOT need to manually instruct subagents about context-mode.

## Output constraints

- Keep responses under 500 words.
- Write artifacts (code, configs, PRDs) to FILES — never return them as inline text. Return only: file path + 1-line description.
- When indexing content, use descriptive source labels so others can `ctx_search(source: "label")` later.

## ctx commands

| Command | Action |
|---------|--------|
| `ctx stats` | Call the `ctx_stats` MCP tool and display the full output verbatim |
| `ctx doctor` | Call the `ctx_doctor` MCP tool, run the returned shell command, display as checklist |
| `ctx upgrade` | Call the `ctx_upgrade` MCP tool, run the returned shell command, display as checklist |
