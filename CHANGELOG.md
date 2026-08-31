# SeqVerify — Changelog

Fork base: `seqverify-mem2/` (local fork of https://github.com/mpiersonsmela/SeqVerify)  
Fork location: `Claude/`  
Date: 2026-06-05  

---

## Bug Fixes

### Java version incompatibility (SnpEff/SnpSift)
**Problem:** The `seqverify` conda environment ships OpenJDK 11, but SnpEff 5.3.0a was compiled against Java 21 (`class file version 65.0 > 55.0`). This caused `UnsupportedClassVersionError` every time SnpEff or SnpSift was invoked, silently producing empty annotation output.

**Fix (`seqverify`):** Replaced `os.system("snpEff ...")` and `os.system("SnpSift ...")` with helper functions `run_snpeff()` / `run_snpsift()` that invoke the JARs directly using the system Java 21 binary (`/usr/lib/jvm/java-21-openjdk-amd64/bin/java -jar`). The bundled SnpEff config (`snpeff-5.3.0a-0/snpEff.config`) is passed with `-config` so that current database names are recognised.

### ClinVar benign filtering
**Problem:** `seqver_lofFinder.py` (and the new `seqver_pathogenicity.py`) replaced underscores with spaces in CLNSIG values before comparing against the `_BENIGN` / `_PATHOGENIC` sets, which store the canonical underscore-delimited ClinVar strings (e.g. `Benign/Likely_benign`). After replacement the strings no longer matched, so ~244 benign variants per run leaked into the TIER2_MODERATE output.

**Fix (`seqver_pathogenicity.py`):** Split CLNSIG on `|` without first replacing underscores. A separate `clnsig_display` variable (underscores → spaces) is used only for the human-readable TSV column.

### Tiering priority (ClinVar Benign vs HIGH impact)
**Problem:** The tiering logic checked `impact == "HIGH"` before `is_benign`, so stop-gained and frameshift variants with a ClinVar Benign classification were reported as TIER1_HIGH_IMPACT instead of being excluded. ClinVar expert curation supersedes SnpEff functional prediction for characterised variants.

**Fix (`seqver_pathogenicity.py`):** Moved the `is_benign → exclude` check to run immediately after the `is_pathogenic` check, before any impact-level tier.

### `commandHandler` called outside pipeline stages
**Problem:** The original script called `commandHandler(..., return_commands_only=False)` unconditionally at the top of `main()`, before any stage guard. This wrote genome files to `temp_folder` even when resuming from late stages (e.g. `start=snp_filtering`), failing if `temp_folder` didn't exist or was read-only.

**Fix (`seqverify`):** Top-level call now uses `return_commands_only=True` (reads the command file, computes coordinate offsets, no I/O). Genome file creation is moved inside the `beginning` stage where it belongs.

### Blank CNVpytor output from stale GC file
**Problem:** `run_cnvpytor_analysis` skipped GC-content generation whenever a GC `.pytor` file already existed (`if not os.path.exists(gc_file_path)`). An empty or truncated GC file left behind by an interrupted run passed this check, so generation was skipped. CNVpytor's `calculate_histograms` only builds histograms for chromosomes present in the GC data, so an empty GC file produced **zero histograms, zero CNV calls, and a blank Manhattan plot — with no error raised**.

**Fix (`seqverify`):** The GC step now validates an existing file with `io_gc.gc_chromosomes() > 0` before reusing it; stale/empty/unreadable files are deleted and regenerated. GC generation also verifies it produced at least one chromosome before the file is trusted.

### Mitochondrial chromosome named `chrMT` not recognized
**Problem:** The chromosome-type classifier in `run_cnvpytor_analysis` matched `chrM`/`MT`/`M`/`mitochondrion` but not `chrMT` (used by some assemblies, e.g. the cynomolgus T2T-MFA8v1.1 genome), so the mitochondrion was mis-typed as an autosome.

**Fix (`seqverify`):** Added `chrMT` to the mitochondrial name list.

---

## New Features

### 1. Liftover-based variant calling (`seqver_liftover.py`)

**Replaces:** Re-alignment of reads to GRCh38 for variant calling.

**How it works:**
1. Variants are called with `bcftools mpileup | bcftools call` directly against the **already-aligned edited-CHM13 BAM** — eliminating a 2–8 hour re-alignment step.
2. `adjust_vcf_for_edits()` remaps variant coordinates from edited-CHM13 space back to unedited-CHM13 space. Each edit's length change is accumulated as a cumulative shift for downstream positions. Variants that fall inside an inserted sequence (no unedited-genome equivalent) are written to a separate `_edit_region.vcf` rather than discarded silently.
3. `liftover_vcf()` calls **CrossMap** (`CrossMap vcf`) with the UCSC hs1→hg38 chain file to lift unedited-CHM13 coordinates to GRCh38. ~95–96% of variants lift successfully; those that fail (CHM13-specific regions) are written to `_unmapped.vcf`.
4. SnpEff annotation then proceeds on GRCh38-coordinate VCF as before.

**New dependency:** CrossMap (`pip install CrossMap`). Install via `setup.sh`.

**New config keys:**
```ini
[VARIANT]
variant_calling=["seqverify_defaults/clinvar.vcf"]   ; ClinVar path only — no GRCh38 FASTA
chain_file=seqverify_defaults/chm13v2.0ToHg38.over.chain.gz
grch38_fasta=seqverify_defaults/Homo_sapiens.GRCh38.dna.primary_assembly.fa  ; optional, for CrossMap allele validation
variant_db=GRCh38.mane.1.2.refseq   ; or GRCh38.p14 etc. — use SnpEff 5.3.x names
```

Note: `variant_calling` is now a **1-element list** (ClinVar VCF only). The GRCh38 FASTA is no longer used for alignment.

**Wild-type path:** If no `--targeted` edit file is provided, `adjust_vcf_for_edits` is a no-op copy and CrossMap runs directly.

### 2. Clinical-grade variant interpretation (`seqver_pathogenicity.py`)

**Replaces:** `seqver_lofFinder.mutation_logger()`.

**Tiering logic (priority order):**
| Tier | Criterion |
|------|-----------|
| `NEAR_EDIT` | Any variant within `variant_window_size` bp of an edit site — always reported regardless of severity |
| `TIER1_CLINVAR_PLP` | ClinVar CLNSIG = Pathogenic / Likely_pathogenic |
| *(excluded)* | ClinVar CLNSIG = Benign / Likely_benign — removed regardless of impact |
| `TIER1_HIGH_IMPACT` | SnpEff impact = HIGH (stop_gained, frameshift, splice_acceptor/donor, start_lost, transcript_ablation, etc.) with no Benign ClinVar entry |
| `TIER2_MODERATE` | SnpEff impact = MODERATE (missense, etc.) with no Benign ClinVar entry |
| *(excluded)* | MODIFIER or LOW impact without P/LP ClinVar support |

**Additional output columns vs original:**
- `CLNSIG` — human-readable ClinVar significance
- `CLNDN` — ClinVar disease name(s)
- `REVIEW_STATUS` — ClinVar review status (confidence proxy)
- `NEAR_EDIT` — True/False
- `DIST_TO_EDIT` — bp distance to nearest edit site
- `TIER` — tiering label above

### 3. CNV quality filtering (`seqver_cnv_filter.py`)

Applied to CNVpytor output in the `snp_filtering` stage (or on demand). Filters CNV calls by:
- `REPEAT_MAPQ` — Q0 > 0.5 (majority of reads have MAPQ=0 → repetitive/centromeric)
- `HIGH_N` — pN > 0.1 (>10% N bases → assembly gap)
- `NOT_SIG` — all three main p-values > 0.05
- `SMALL` — size < 3 × bin_size (below reliable detection threshold)
- `REPEAT_BED` — optional overlap with a user-supplied BED of blacklist/repetitive regions

Outputs: annotated TSV with `FILTER` column + a `_pass.tsv` containing only `PASS` calls.

### 4. LOH detection around edit sites (`seqver_loh.py`)

For each targeted edit site, examines a configurable window (default 1 Mb) in the lifted-over GRCh38 VCF for evidence of loss of heterozygosity:
1. Collects heterozygous calls (`GT=0/1`) with allele depth (`AD`) in the window.
2. Computes allele balance (AB = alt / total) at each het site.
3. Applies a two-sided binomial test (expected AB = 0.5) on pooled read counts.
4. Reports `LOH_CANDIDATE` if p < 0.01 **and** mean AB < 0.3 or > 0.7; `BALANCED` otherwise; `INSUFFICIENT_DATA` if fewer than 5 het SNPs in window.

Output TSV: `EDIT_SITE, CHR, WINDOW_START, WINDOW_END, N_HET_SNPS, MEAN_AB, SUM_ALT, SUM_TOTAL, BINOMIAL_PVAL, LOH_STATUS`.

Wild-type runs (no edit file): produces an empty-body TSV with header only.

**New config key:**
```ini
[VARIANT]
loh_window=1000000   ; half-window in bp
```

### 5. Work directory support (`seqverify`)

New `work_dir` config key / `--work_dir` argument. When set, the pipeline changes to that directory before creating output folders, so all outputs (folders, VCF files, TSVs) are created under `work_dir`. Input paths (reads, genome, databases) are converted to absolute paths before the `chdir`. Useful for separating output from the working genome data directory.

```ini
[MAIN]
work_dir=/path/to/output_location
```

### 6. SnpEff database name updated

The original pipeline hardcoded `GRCh38.105` (an Ensembl release number valid for SnpEff 5.1). SnpEff 5.3.0a does not include this database. The default is now `GRCh38.mane.1.2.refseq` (MANE Select + Plus Clinical transcripts, RefSeq IDs — the most clinically relevant choice for 2025+). The database name is configurable via `variant_db=` in `[VARIANT]`.

---

## New Files

| File | Purpose |
|------|---------|
| `seqver_liftover.py` | CHM13→GRCh38 coordinate adjustment + CrossMap wrapper |
| `seqver_pathogenicity.py` | Clinical variant tiering (replaces `seqver_lofFinder.py`) |
| `seqver_cnv_filter.py` | CNVpytor call quality filtering |
| `seqver_loh.py` | LOH detection at edit sites |
| `setup.sh` | Install CrossMap dependency |
| `CHANGELOG.md` | This file |
| `benchmark/create_benchmark.sh` | Synthetic yeast benchmark for fast pipeline testing |
| `benchmark/yeast.config.template` | Config template for yeast benchmark |

## Unchanged Files (copied verbatim from base)

`seqver_functions.py`, `seqver_plots.py`, `seqver_gtf.py`, `seqver_genomeupdate.py`, `seqver_lofFinder.py` (kept for backward compatibility; superseded by `seqver_pathogenicity.py`).
