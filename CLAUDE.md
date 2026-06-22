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

**Entrypoint:** `cli.py` is a `click` group; the `run` command orchestrates the pipeline. Required options: `--input` (FASTQ), `--output` (dir), `--db` (Tier-1 DB), `--host-ref`, `--specimen`. Notable extras: `--read-type short|long`, `--amr-db` (ABRicate DB, default `card`), `--no-pdf`, `--sketch-threshold`, `--n-bootstrap`.

**Pipeline flow:**
```
QC (fastp/Chopper) → Host removal (BWA-MEM2/minimap2) → Sketch screen (sourmash)
  → Targeted alignment (minimap2) → EM abundance → AMR screen (ABRicate)
  → Report (JSON + TSV + PDF + HTML)
```

**Key data contracts between stages:**
- `PipelineConfig` (`config.py`) — single config dataclass passed to every stage; holds paths, `ReadType` (short/long), `SpecimenType` (blood/csf/bal/tissue), and tuning params
- `EMResult` (`em.py`) — `abundances` array (sums to 1) + `n_reads`; produced by `em_abundance()` and consumed by `write_report()`
- `AlignResult` (`align.py`) — `alignment_matrix` (reads × organisms binary array) + `organism_names`; bridges alignment → EM
- `ReportEntry` (`report.py`) — one row per organism with lazy `grade` property; contaminant flag set externally by `flag_contaminants()`
- `AMRHit` (`amr.py`) — `run_amr_screen()` converts reads to FASTA and runs ABRicate, mapping resistance genes back to `organism_names`; results are attached to the report
- Report rendering: `report.py` writes JSON + TSV; `pdf_report.py` (`write_pdf_report`) and `html_report.py` (`write_html_report`) render the same `ReportEntry` rows into clinical PDF/HTML

**Evidence grading** (`report.py` + `contaminants.py`):
- Grade A/B/C/X is a property on `ReportEntry`, not a stored field
- Grading uses specimen-specific minimum read thresholds (`_MIN_READS`) + CI width ≤ 0.15 for Grade A
- `contaminants.py` holds `CONTAMINANT_PRIORS` — a per-`SpecimenType` list of common contaminant species; `flag_contaminants()` sets `contaminant_risk=True` via case-insensitive substring match, which demotes Grade A → B/C

**Read-type branching:**
- Short reads: fastp for QC, BWA-MEM2 for host removal
- Long reads: Chopper for QC, minimap2 for host removal
- Both paths produce identical intermediate files; downstream stages are read-type agnostic

**External tool dependencies:** fastp, Chopper, BWA-MEM2, minimap2, samtools, sourmash, abricate — must be on `PATH` (install via conda bioconda channel)

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
