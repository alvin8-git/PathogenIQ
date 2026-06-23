# Open-world detection for airborne pathogens (incl. novel)

The targeted arm (sketch → align → EM vs the ~110-genome tier-1 DB) is
**closed-world**: blind to any pathogen not in the DB, including novel ones. Air
surveillance must catch (a) known pathogens off-DB and (b) genuinely novel ones,
and air is disproportionately **viral**. This is the requirements map and status.

| # | Requirement | Role | Status |
|---|-------------|------|--------|
| R1 | **Novelty trigger** (`novelty.py`) — Kraken2 vs a *broad* DB; report unclassified "dark-matter" fraction | Cheap gate: "is there something with no reference here?" | **built** (`--novelty`) |
| R2 | **Viral arm** (`viral.py`) — geNomad (ID + ICTV taxonomy) → CheckV (completeness) on contigs | Air pathogens are viral-heavy; bacterial/MAG arms can't see them | **built** (`--viral`) |
| R3 | **Novel-bacterial recovery** — assemble → MetaBAT2 → GTDB-Tk placement | Reference-free recovery of unknown bacteria | built (`--assemble`) |
| R4 | **Pathogenicity triage** — contig-level VFDB/CARD + phylo-proximity-to-pathogen | Discriminator: novel *pathogen* vs novel benign environmental microbe | **next** |
| R5 | **Open-world grading** — evidence tier for reference-free hits (breadth, completeness, hallmark/marker count, phylo-confidence) | Consistent A/B/C/X-equivalent for assembled/novel/viral hits | **next** |

## Design notes (R1, R2 — built 2026-06-23)

- **Broad DB, not tier-1.** Novelty classifies against `databases/kraken2`
  (Standard); the tier-1 kraken DB would mark all non-target reads unclassified
  and falsely inflate novelty.
- **Runs before the no-hits early return.** Novelty + the open-world arms now
  execute even when the targeted screen finds zero candidates — a viral-dominant
  or wholly-novel air sample no longer exits empty. The CLI writes a report with
  empty `findings` plus the `novelty` / `mags` / `viral` blocks.
- **Non-blocking.** Every arm degrades to a skip when its tool or DB is absent
  (kraken2 + Standard DB; geNomad + `$GENOMAD_DB`; CheckV + `$CHECKV_DB`), matching
  the AMR/assembly pattern. `environment.yml` pins genomad + checkv; kraken2 was
  already present.
- **Source caveat (viral).** DNA mNGS captures DNA viruses + integrated proviruses,
  NOT RNA viruses (influenza, SARS-CoV-2, RSV). The arm is source-agnostic — feed
  it metatranscriptome contigs and it classifies RNA viruses too; the capture is a
  wet-lab (DNA vs RNA extraction) decision, recorded in the wet-lab section of `todo.md`.

## Validating the viral arm (R2)

Our air data can't validate it (DNA-seq → no RNA viruses; Zymo D6300 has no viral
members; filters carry unlabelled phage only). `scripts/12_viral_insilico_spikein.py`
generates gold truth instead: simulate reads (wgsim) from known viral genomes
(T4 + Lambda phage, SARS-CoV-2) → optionally mix into a real aircraft-filter
background → megahit → `run_viral_stage` → score geNomad recall (correct ICTV
lineage) + CheckV completeness vs truth. Scoring logic is offline-verifiable
(`--selfcheck`); the full run needs wgsim + megahit + geNomad + CheckV (+ DBs).
This validates the *bioinformatic* arm; end-to-end air-virus capability still needs
real RNA-seq air virome (a wet-lab decision).

### Result — first scored run (2026-06-23): RECALL 100% (3/3)

30× wgsim reads from each truth genome mixed into the SRR32514319 aircraft-filter
background (~6.5 M reads), assembled (megahit → 56,442 contigs), geNomad + CheckV.
geNomad called **190 viral contigs** (most are real background phage — bonus); all
three spiked genomes recovered with correct ICTV lineage:

| Genome | Expected | Recovered | CheckV completeness | geNomad lineage |
|--------|----------|-----------|---------------------|-----------------|
| Escherichia phage T4 | Straboviridae | ✅ | 21.6% | …Caudoviricetes;Straboviridae |
| Escherichia phage Lambda | Caudoviricetes | ✅ | 76.2% | …Caudoviricetes |
| SARS-CoV-2 | Coronaviridae | ✅ | 70.2% | …Nidovirales;Coronaviridae |

Notes: detection (recall) is the headline — all 3 found, lineage correct, including
the respiratory RNA virus classified from its genome amid real air-metagenome noise.
T4's low completeness (21.6%) reflects megahit fragmenting its 169 kb genome at 30×,
not a miss. This run also surfaced+fixed a real bug: `run_megahit` failed silently
when `assembly/` didn't pre-exist (megahit's `os.mkdir` is single-level), which had
made `--assemble` and `--viral` no-op on every run — caught only because the harness
expected ≥3 contigs and got 0.

**Tooling:** geNomad + CheckV live in a dedicated `genomad` conda env (numpy<2),
CLIs symlinked into the pathogeniq env (`scripts/13_setup_viral_env.sh`). DBs at
`databases/genomad_db` (v1.9) + `databases/checkv_db` (v1.5).

## JSON additions

- `novelty`: `{total_reads, classified_reads, unclassified_reads, unclassified_fraction, n_species, top_taxa[], flagged}`
- `viral`: `[{contig_id, length, taxonomy, topology, virus_score, n_hallmarks, completeness, checkv_quality}]`

## Next: R4 pathogenicity triage

The earlier air goal — "detect pathogens, not environmental organisms" — is R4.
A novel contig/MAG matters as a *pathogen* only if it carries virulence (VFDB) /
AMR (CARD) / toxin genes **or** sits phylogenetically adjacent to a known pathogen
(e.g. a novel *Bacillus* near *B. anthracis*). Plan: run VFDB/CARD on contigs (not
just reads) + a phylo-proximity score off the GTDB-Tk placement, feeding R5's
open-world grade.
