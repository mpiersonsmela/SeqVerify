# SeqVerify
SeqVerify is a Python-based command line tool for analysis of whole genome sequencing data for gene-editing verification. It performs insertion site detection, copy number variation (CNV) analysis through CNVPytor, bacterial contamination detection through KRAKEN2 and BRACKEN, and variant calling and filtering aided by SnpEff and SnpSift.

## Install

### Dependencies
* BCFtools
* BLAST
* BRACKEN
* BWA ≥0.7 (or BWA-MEM2 with `--use_mem2`)
* CNVPytor ≥1.3
* CrossMap
* HTSLIB
* IDNA
* IGV ≥2.13.2
* IGVReports
* Kraken2 ≥2.0
* Matplotlib ≥2.2
* Numpy
* **OpenJDK 21** (required for SnpEff 5.3+; older Java versions will fail)
* Pysam (non-HDR allele classifier)
* Python ≥3.10
* SAMtools ≥1.14
* SciPy
* SNPEff ≥5.3
* SNPSift ≥5.3
* SPAdes (de-novo assembly of non-HDR reads)

### Install through Bioconda

```
conda install -c bioconda seqverify
```

## Usage

```
seqverify --output output_name --reads_1 sample_1.fastq --reads_2 sample_2.fastq \
          --genome genome.fa --untargeted transgenes.fa --targeted commands.txt \
          --database db8gb
```

A configuration file is recommended for complex runs — see **Config File** below.

### Anatomy of a SeqVerify call

SeqVerify has the following standard arguments:
* `--output` (String) Identifying name for all output files and folders. E.g. `--output Sample1` creates `Sample1_seqverify/`.
* `--reads_1` / `--reads_2` (String/Path) Paired-end FASTQ (or gzipped FASTQ) files.
* `--genome` (String/Path) Reference genome FASTA. Defaults to CHM13v2.0 (see `--download_defaults`).
* `--untargeted` (String/Path) FASTA file(s) containing sequences of markers to detect (transgenes, unwanted plasmids, etc.). Space-separated; compatible with `--targeted`.
* `--targeted` (String/Path) Command file specifying known insertion sites. See **Command text file** below.
* `--gtf` (String/Path) GTF/GFF3 annotation file, updated with edits from `--targeted`.

#### Performance
* `--threads` (Integer) CPU threads. Default: 1.
* `--max_mem` (String) Maximum memory for genome indexing (e.g. `16G`). Default: `16G`.
* `--use_mem2` Use BWA-MEM2 instead of BWA for alignment (faster).
* `--start` (String) Resume from a specific pipeline stage:
  `all` (default), `beginning`, `align`, `markers`, `readout`, `nonhdr`, `cnv`, `plots`, `kraken`, `variant`, `snp_filtering`.
* `--deduplicate` / `--no-deduplicate` Mark PCR/optical duplicates (`samtools markdup`) before CNV and variant calling so they do not inflate read depth or allele fractions. Enabled by default. (Config: `[INSERTION] deduplicate`.)
* `--analyze_nonhdr` / `--no-analyze-nonhdr` Run the non-HDR allele analysis at targeted edit sites (see below). Enabled by default; skipped automatically when there is no `--targeted` edits file. (Config: `[INSERTION] analyze_nonhdr`.)

#### KRAKEN2
* `--kraken` Enable KRAKEN2 microbial contamination analysis.
* `--database` (String/Path) Path to a KRAKEN2 database.

#### Insertion Site Detection
* `--granularity` (Integer) Bin size for merging nearby insertion calls (bp). Default: 500.
* `--min_matches` (Integer) Minimum reads per insertion site to appear in the readout. Default: 1.
* `--mitochondrial` Detect mitochondrial DNA insertions (chrM).

#### CNVPytor
* `--bin_size` (Integer) Bin size for CNV detection (bp). Default: 100,000.
* `--manual_plots` Use matplotlib instead of IGVReports for insertion site coverage plots.

#### Variant Calling
Variants are called from the CHM13-aligned BAM (reusing the main alignment — no second alignment to GRCh38 needed), then lifted over to GRCh38 coordinates via CrossMap before annotation.

* `--variant_calling` Path to a ClinVar VCF for annotation. Takes **one** argument (unlike older versions):
  ```
  --variant_calling seqverify_defaults/clinvar.vcf
  ```
* `--chain_file` Path to a CHM13v2.0 → hg38 CrossMap chain file (downloaded by `--download_defaults`).
* `--grch38_fasta` GRCh38 FASTA for CrossMap allele validation (optional but recommended).
* `--variant_db` SnpEff database name. Default: `GRCh38.mane.1.2.refseq` (SnpEff 5.3+).
* `--variant_intensity` Minimum SnpEff impact level to report (`MODIFIER`, `LOW`, `MODERATE`, `HIGH`). Default: `MODERATE`.
* `--variant_window_size` (Integer) Distance from an edit site within which any variant is flagged as `NEAR_EDIT` regardless of severity. Default: 10,000 bp.
* `--loh_window` (Integer) Half-window size around each edit site for LOH detection (bp). Default: 1,000,000.
* `--min_quality` (Integer) Minimum bcftools quality score to include a variant. Default: 3.

#### Variant output
`seqverify_output_variants.tsv` reports variants with the following tiering:

| Tier | Criteria |
|------|----------|
| `TIER1_CLINVAR_PLP` | ClinVar Pathogenic / Likely_pathogenic |
| `TIER1_HIGH_IMPACT` | SnpEff HIGH impact (stop_gained, frameshift, splice donor/acceptor, start_lost) |
| `TIER2_MODERATE` | SnpEff MODERATE impact, ClinVar not Benign |
| `NEAR_EDIT` | Within `--variant_window_size` of a targeted edit site |
| *(excluded)* | ClinVar Benign/Likely_benign; MODIFIER impact without ClinVar support |

`seqverify_output_loh.tsv` reports allele-balance statistics for each edit site window. Sites with a binomial test p < 0.01 and mean allele balance outside [0.3, 0.7] are flagged `LOH_CANDIDATE`.

`copy_number/calls.{bin_size}_filtered.tsv` adds a `FILTER` column to CNVpytor calls (flags: `REPEAT_MAPQ`, `HIGH_N`, `NOT_SIG`, `SMALL`). `calls_{bin_size}_filtered_pass.tsv` contains only `PASS` calls.

#### Other
* `--keep_temp` Keep temporary files. Off by default (temp files can exceed 100 GB).
* `--download_defaults` Download default reference files: CHM13v2.0, GRCh38, Kraken PlusPFP 8GB, ClinVar, SnpEff config, and the CHM13v2.0→hg38 CrossMap chain file.

## Output

SeqVerify creates two directories:
* `output_seqverify/` — final results
* `output_seqverify_temp/` — intermediate files (deleted unless `--keep_temp` is set)

Key output files in `output_seqverify/insertion/`:
* `seqverify_readout.txt` — insertion site readout (chromosome, position, marker, read counts, confidence score)
* `igv_viewer.html` — IGV-based coverage viewer for insertion sites

In `output_seqverify/copy_number/`:
* `calls.{bin_size}.tsv` — raw CNVpytor calls
* `calls.{bin_size}_filtered.tsv` / `_filtered_pass.tsv` — quality-filtered CNV calls

In `output_seqverify/variant_calling/`:
* `seqverify_output.ann.vcf` — SnpEff + ClinVar annotated VCF (GRCh38 coordinates)
* `seqverify_output_variants.tsv` — clinically tiered variant table
* `seqverify_output_loh.tsv` — LOH analysis at edit sites

## Config File

A template config file (`seqverify.config`) is bundled with SeqVerify. Use it with `--config`:

```
seqverify --config seqverify.config
```

Do not modify the section headers (text in square brackets). File names in `untargeted` must not contain spaces.

Example `[VARIANT]` section:
```ini
[VARIANT]
variant_calling=["seqverify_defaults/clinvar.vcf"]
chain_file=seqverify_defaults/chm13v2.0ToHg38.over.chain.gz
grch38_fasta=seqverify_defaults/Homo_sapiens.GRCh38.dna.primary_assembly.fa
variant_db=GRCh38.mane.1.2.refseq
variant_intensity=MODERATE
variant_window_size=10000
loh_window=1000000
```

## Command text file

The `--targeted` command file specifies modifications to the reference genome at known edit sites. Format (tab-separated):

```
CHR:START-END	SEQUENCE
```

* `CHR:START-END` — chromosome and 0-based Python coordinates of the region to replace
* `SEQUENCE` — replacement sequence

Example: `chr14:54643109-54643153	CTAGATATCGGCGCGCC...` replaces 45 bp at chr14:54,643,109–54,643,153 with the provided sequence.

Multiple commands per file are supported. SeqVerify automatically resolves coordinate shifts when multiple edits affect the same chromosome.

## Spurious Filtering

`--spurious_filtering_threshold` (default: 0.00001) controls filtering of extremely high-coverage insertion sites. Set to `0` to disable.

## Non-HDR allele analysis

For each `--targeted` edit site, SeqVerify reconstructs and annotates the allele that did
**not** receive the HDR cassette (enabled by default; `--no-analyze-nonhdr` to skip, requires
SPAdes). Read pairs spanning the locus are realigned to a wild-type ±10 kb mini-reference —
HDR reads soft-clip at the cassette junction while non-HDR reads span the cut — the non-HDR
reads are separated and de-novo assembled with SPAdes, and the result is annotated for CDS
impact and distance to the nearest splice site. Outputs land in
`output_seqverify/nonHDR_allele/`:
* `SUMMARY.txt` — per-site scar description, gene/CDS/splice annotation, and a zygosity verdict
* `<site>_nonHDR_igv.html` — standalone IGV report (reads realigned to WT + assembled contig + gene model)
* `<site>_wtaln.bam`, `<site>_spades/contigs.fasta`, `<site>_ctg_vs_wt.bam` — intermediates

## Helper: generate an edits file from an HDR donor plasmid

`scripts/make_edits_from_plasmid.py` builds a `--targeted` command file from a GenBank HDR
donor plasmid by BLASTing its homology arms against the reference and emitting the cargo
between them. Rox-flanked selection cassettes are collapsed to a single Rox scar by default
(`--keep-rox` disables this). Uses NCBI BLAST+ (already a dependency); no Biopython.

```
python scripts/make_edits_from_plasmid.py --plasmid donor.gb --ref chm13v2.0.fa --out edits.txt
```

## Tests

```
python tests/test_seqverify.py     # or: pytest tests/
```

## Frequently Asked Questions

### Why do I get a Java version error with SnpEff?

SnpEff 5.3+ requires **Java 21**. The conda environment now pins `openjdk=21`. If you see `UnsupportedClassVersionError: class file version 65.0`, update your environment:
```
conda env update -f seqverify-env.yml
```

### Why does --variant_calling now take only one argument?

Variant calling was redesigned to reuse the CHM13 alignment and lift variants to GRCh38 via CrossMap, eliminating the need for a second alignment to GRCh38. Pass only the ClinVar VCF path. The chain file is specified via `--chain_file`.

### Why am I getting a BCFTools error?

If you see `bcftools: error while loading shared libraries: libgsl.so.25`, try forcing bioconda in your `.condarc`:
```yaml
channel_priority: strict
channels:
  - bioconda
  - conda-forge
  - defaults
```

### Can multiple FASTQ files be used?

No; concatenate them first with `cat`.

### Are insertions of mitochondrial DNA a sign of contamination?

Not necessarily — NUMTs (nuclear mitochondrial DNA segments) are common in the human genome.

## Citation

If you use SeqVerify, please cite:

> Smela MP, Pepe V, Lubbe S, Kiskinis E, Church GM. SeqVerify: An accessible analysis tool for cell line genomic integrity, contamination, and gene editing outcomes. *Stem Cell Reports*. 2024;19(10):1505–1515. doi:[10.1016/j.stemcr.2024.08.004](https://doi.org/10.1016/j.stemcr.2024.08.004)

Citation metadata is also provided in [`CITATION.cff`](CITATION.cff); GitHub's "Cite this repository" (on the right side of root) button generates APA and BibTeX from it.
