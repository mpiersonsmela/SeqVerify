"""
seqver_loh.py
Detect loss of heterozygosity (LOH) in windows around edit sites.

Strategy:
  1. Parse SnpEff+ClinVar VCF for heterozygous SNPs (GT 0/1 or 1/0).
  2. For each edit site, collect het SNPs within a user-defined window.
  3. Compute allele balance (AB = alt_depth / total_depth) at each het site.
  4. Use a two-sided binomial test (expected AB = 0.5) on the pooled read counts.
  5. Flag windows with p < 0.01 and mean AB far from 0.5 as LOH_CANDIDATE.

Output TSV columns:
  EDIT_SITE  CHR  WINDOW_START  WINDOW_END
  N_HET_SNPS  MEAN_AB  SUM_ALT  SUM_TOTAL  BINOMIAL_PVAL  LOH_STATUS

LOH_STATUS values:
  LOH_CANDIDATE        p < pval_threshold and mean AB outside ab_range
  BALANCED             het SNPs present, no significant skew
  INSUFFICIENT_DATA    fewer than min_het_snps het SNPs in window
  NO_EDIT_SITES        wild-type run — output file written with empty body
"""

import os


def _parse_het_variants(vcf_file):
    """
    Parse VCF and return list of (chrom, pos, alt_depth, ref_depth) for
    heterozygous calls that have AD (allele depth) information.
    """
    het_variants = []
    with open(vcf_file) as f:
        fmt_col = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.lstrip("#").strip().split("\t")
                try:
                    fmt_col = cols.index("FORMAT")
                except ValueError:
                    fmt_col = 8  # standard position
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            fmt_keys = fields[8].split(":")
            samp_vals = fields[9].split(":")

            if len(fmt_keys) != len(samp_vals):
                continue

            fmt = dict(zip(fmt_keys, samp_vals))
            gt = fmt.get("GT", "./.")
            alleles = gt.replace("|", "/").split("/")

            # Heterozygous: exactly two distinct alleles, neither is '.'
            if len(set(alleles)) != 2 or "." in alleles:
                continue

            ad_str = fmt.get("AD", "")
            if not ad_str or "," not in ad_str:
                continue

            try:
                depths = [int(x) for x in ad_str.split(",")]
                ref_d = depths[0]
                alt_d = sum(depths[1:])  # sum alt allele depths if multi-allelic
            except ValueError:
                continue

            total = ref_d + alt_d
            if total < 4:  # skip very low coverage sites
                continue

            het_variants.append((chrom, pos, alt_d, ref_d))

    return het_variants


def _build_edit_sites(edit_commands, command_file):
    """Return list of (chrom, center, label) for each edit site."""
    sites = []
    if edit_commands:
        for cmd in edit_commands:
            center = (cmd.coords[0] + cmd.coords[1]) // 2
            label = f"{cmd.chr}:{cmd.coords[0]}-{cmd.coords[1]}"
            sites.append((cmd.chr, center, label))
    elif command_file and os.path.exists(command_file):
        with open(command_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")[0].split(":")
                chrom = parts[0]
                coords = [int(x) for x in parts[1].split("-")]
                center = (coords[0] + coords[1]) // 2
                label = f"{chrom}:{coords[0]}-{coords[1]}"
                sites.append((chrom, center, label))
    return sites


def detect_loh_around_edits(vcf_file, edit_commands, output_file,
                            command_file=None,
                            window=1_000_000,
                            min_het_snps=5,
                            pval_threshold=0.01,
                            ab_range=(0.3, 0.7)):
    """
    Detect LOH around each edit site.

    vcf_file:      path to final annotated VCF (GRCh38 coordinates)
    edit_commands: list of editArgument objects, or None for wild-type
    output_file:   path for output TSV
    command_file:  fallback if edit_commands is None
    window:        half-window size in bp around each edit center
    min_het_snps:  minimum het SNPs required to make a call
    pval_threshold: binomial test significance threshold
    ab_range:      (low, high) — mean AB outside this range is LOH-consistent
    """
    try:
        from scipy.stats import binomtest
    except ImportError:
        # scipy < 1.7 uses binom_test
        from scipy.stats import binom_test as _bt
        def binomtest(k, n, p):
            class _R:
                pvalue = _bt(k, n, p)
            return _R()

    edit_sites = _build_edit_sites(edit_commands, command_file)

    header = (
        "EDIT_SITE\tCHR\tWINDOW_START\tWINDOW_END\t"
        "N_HET_SNPS\tMEAN_AB\tSUM_ALT\tSUM_TOTAL\t"
        "BINOMIAL_PVAL\tLOH_STATUS\n"
    )

    with open(output_file, "w") as out:
        out.write(header)

        if not edit_sites:
            # Wild-type genome: no edit sites to check
            print("LOH detection: no edit sites — writing empty output")
            return

        if not os.path.exists(vcf_file):
            print(f"LOH detection: VCF not found: {vcf_file}")
            return

        het_vars = _parse_het_variants(vcf_file)
        print(f"LOH detection: {len(het_vars)} heterozygous SNPs parsed from VCF")

        # Index het variants by chromosome for fast lookup
        by_chrom = {}
        for chrom, pos, alt_d, ref_d in het_vars:
            by_chrom.setdefault(chrom, []).append((pos, alt_d, ref_d))

        for site_label, chrom, center in [(s[2], s[0], s[1]) for s in edit_sites]:
            win_start = max(0, center - window)
            win_end = center + window

            # Collect het SNPs in window
            snps = [(pos, ad, rd) for pos, ad, rd in by_chrom.get(chrom, [])
                    if win_start <= pos <= win_end]

            n = len(snps)

            if n < min_het_snps:
                status = "INSUFFICIENT_DATA"
                mean_ab = "."
                sum_alt = sum_total = 0
                pval = "."
            else:
                sum_alt = sum(ad for _, ad, rd in snps)
                sum_total = sum(ad + rd for _, ad, rd in snps)
                mean_ab = sum_alt / sum_total if sum_total > 0 else 0.5

                result = binomtest(sum_alt, sum_total, 0.5)
                pval = result.pvalue

                loh_consistent = mean_ab < ab_range[0] or mean_ab > ab_range[1]
                if pval < pval_threshold and loh_consistent:
                    status = "LOH_CANDIDATE"
                else:
                    status = "BALANCED"

                mean_ab = f"{mean_ab:.4f}"
                pval = f"{pval:.2e}"

            row = [
                site_label,
                chrom,
                str(win_start),
                str(win_end),
                str(n),
                str(mean_ab),
                str(sum_alt),
                str(sum_total),
                str(pval),
                status,
            ]
            out.write("\t".join(row) + "\n")
            print(f"  {site_label}: {n} het SNPs, status={status}")
