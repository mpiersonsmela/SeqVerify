"""
seqver_cnv_filter.py
Filter CNVpytor TSV calls by quality metrics and optionally by overlap with
repetitive / low-mappability regions.

CNVpytor call columns (from seqverify's run_cnvpytor_analysis):
  type  chrom  start  end  size  cnv  p_val  p_val_2  p_val_3  p_val_4  Q0  pN  dG

Filter criteria:
  REPEAT_MAPQ  Q0 > 0.5        (>50% zero-MAPQ reads — repetitive region)
  HIGH_N       pN > 0.1        (>10% N bases — assembly gap / low-quality)
  NOT_SIG      all three main p-values > 0.05
  SMALL        size < 3×bin_size
  REPEAT_BED   overlaps provided repetitive-region BED (optional)
  PASS         none of the above

Outputs:
  output_tsv         all calls with added FILTER column
  output_pass_tsv    PASS calls only (same path with _pass suffix)
"""

import os


def _parse_bed(bed_file):
    """Load BED into {chrom: [(start, end), ...]} sorted by start."""
    regions = {}
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.setdefault(chrom, []).append((start, end))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def _overlaps_bed(chrom, start, end, bed_regions):
    """Return True if (chrom, start, end) overlaps any region in bed_regions."""
    if chrom not in bed_regions:
        return False
    for rs, re in bed_regions[chrom]:
        if rs >= end:
            break
        if re > start:
            return True
    return False


def filter_cnv_calls(calls_tsv, output_tsv,
                     repetitive_bed=None, bin_size=100000,
                     q0_threshold=0.5, pn_threshold=0.1,
                     pval_threshold=0.05):
    """
    Apply quality filters to CNVpytor calls TSV and write annotated output.

    calls_tsv:      path to calls.<bin_size>.tsv from run_cnvpytor_analysis
    output_tsv:     path for annotated output (FILTER column appended)
    repetitive_bed: optional BED file of repetitive/blacklist regions
    bin_size:       CNVpytor bin size used (for SMALL filter)
    """
    if not os.path.exists(calls_tsv):
        print(f"CNV filter: {calls_tsv} not found, skipping")
        return

    bed_regions = _parse_bed(repetitive_bed) if repetitive_bed else {}
    pass_path = output_tsv.replace(".tsv", "_pass.tsv")
    if not pass_path.endswith("_pass.tsv"):
        pass_path = output_tsv + "_pass.tsv"

    total = 0
    passed = 0
    filter_counts = {}

    with open(calls_tsv) as fin, \
         open(output_tsv, "w") as fall, \
         open(pass_path, "w") as fpass:

        header = fin.readline()
        col_headers = header.rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(col_headers)}

        annotated_header = header.rstrip("\n") + "\tFILTER\n"
        fall.write(annotated_header)
        fpass.write(annotated_header)

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            total += 1

            reasons = []

            try:
                chrom = fields[col["chrom"]]
                start = int(fields[col["start"]])
                end = int(fields[col["end"]])
                size = int(fields[col["size"]])
                q0 = float(fields[col["Q0"]]) if "Q0" in col else 0.0
                pn = float(fields[col["pN"]]) if "pN" in col else 0.0
                p1 = float(fields[col["p_val"]]) if "p_val" in col else 1.0
                p2 = float(fields[col["p_val_2"]]) if "p_val_2" in col else 1.0
                p3 = float(fields[col["p_val_3"]]) if "p_val_3" in col else 1.0
            except (IndexError, ValueError, KeyError):
                reasons.append("PARSE_ERROR")
                fall.write(line + "\t" + ";".join(reasons) + "\n")
                continue

            if q0 > q0_threshold:
                reasons.append("REPEAT_MAPQ")
            if pn > pn_threshold:
                reasons.append("HIGH_N")
            if p1 > pval_threshold and p2 > pval_threshold and p3 > pval_threshold:
                reasons.append("NOT_SIG")
            if size < 3 * bin_size:
                reasons.append("SMALL")
            if bed_regions and _overlaps_bed(chrom, start, end, bed_regions):
                reasons.append("REPEAT_BED")

            filter_val = ";".join(reasons) if reasons else "PASS"
            for r in (reasons or ["PASS"]):
                filter_counts[r] = filter_counts.get(r, 0) + 1

            annotated = line + "\t" + filter_val + "\n"
            fall.write(annotated)
            if not reasons:
                fpass.write(annotated)
                passed += 1

    print(f"CNV filter: {total} calls → {passed} PASS")
    for reason, count in sorted(filter_counts.items()):
        print(f"  {reason}: {count}")
    print(f"  PASS calls written to: {pass_path}")
