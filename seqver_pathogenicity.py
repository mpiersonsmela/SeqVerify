"""
seqver_pathogenicity.py
Clinical-grade tiered variant filtering replacing seqver_lofFinder.mutation_logger().

Tiering logic:
  TIER1  ClinVar P/LP                                → always report
  TIER1  SnpEff HIGH impact                          → always report
  TIER2  SnpEff MODERATE + ClinVar not B/LB          → report
  NEAR_EDIT  variant within window of edit site       → report regardless of tier
  EXCLUDE  ClinVar B/LB (and not NEAR_EDIT)          → drop
  EXCLUDE  MODIFIER impact + no ClinVar entry        → drop

Output TSV columns:
  CHR:COORD  QUALITY  VARIANT_TYPE  EFFECT  GENE  REFSEQ_ID
  DNA_change  PROTEIN_change  HOMOZYGOUS  LOF
  CLNSIG  CLNDN  REVIEW_STATUS  NEAR_EDIT  TIER
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ClinVar significance categories
_PATHOGENIC = {"Pathogenic", "Likely_pathogenic",
               "Pathogenic/Likely_pathogenic"}
_BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
_UNCERTAIN = {"Uncertain_significance", "Conflicting_interpretations_of_pathogenicity",
              "not_provided", ""}

_HIGH_EFFECTS = {
    "stop_gained", "frameshift_variant", "splice_acceptor_variant",
    "splice_donor_variant", "start_lost", "stop_lost",
    "transcript_ablation", "exon_loss_variant", "gene_fusion",
}


def _parse_clnsig(info_fields):
    """Extract CLNSIG and CLNDN from a dict of INFO key=value pairs."""
    clnsig_raw = info_fields.get("CLNSIG", "")
    clndn = info_fields.get("CLNDN", "").replace("_", " ")
    revstat = info_fields.get("CLNREVSTAT", "").replace("_", " ")
    # Keep underscores when splitting — ClinVar stores them as e.g. "Benign/Likely_benign".
    # The sets _PATHOGENIC, _BENIGN, _UNCERTAIN all use the underscore form.
    clnsig_vals = set(clnsig_raw.split("|"))
    if clnsig_vals & _PATHOGENIC:
        tier_clnsig = "Pathogenic_or_LP"
    elif clnsig_vals & _BENIGN:
        tier_clnsig = "Benign_or_LB"
    elif clnsig_vals - {""} - _UNCERTAIN:
        tier_clnsig = next(iter(clnsig_vals - {""} - _UNCERTAIN))
    else:
        tier_clnsig = "None"
    # Human-readable display for the CLNSIG output column
    clnsig_display = clnsig_raw.replace("_", " ")
    return tier_clnsig, clndn, revstat, clnsig_display


def _parse_info(info_str):
    """Parse VCF INFO field into a dict."""
    d = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        else:
            d[token] = True
    return d


def _build_edit_coords(command_file, edit_commands_list):
    """Return {chrom: [start_pos, ...]} from edit command objects or a file path."""
    coords = {}
    if edit_commands_list:
        for cmd in edit_commands_list:
            coords.setdefault(cmd.chr, []).append(cmd.coords[0])
    elif command_file:
        with open(command_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")[0].split(":")
                chrom = parts[0]
                start = int(parts[1].split("-")[0])
                coords.setdefault(chrom, []).append(start)
    return coords


def pathogenicity_logger(folder, output_file, vcf_file,
                         min_quality, min_intensity="MODERATE",
                         window_size=10000,
                         command_file=None, edit_commands=None):
    """
    Parse an SnpEff+ClinVar-annotated VCF and write a clinically-tiered TSV.

    edit_commands: list of editArgument objects (preferred) or None.
    command_file:  path to targeted edits file (fallback if edit_commands is None).
    """
    intensities = ["MODIFIER", "LOW", "MODERATE", "HIGH"]
    min_idx = intensities.index(min_intensity)

    edit_coords = _build_edit_coords(command_file, edit_commands)
    qual_scores = []

    header = (
        "CHR:COORD\tQUALITY\tVARIANT_TYPE\tEFFECT\tGENE\tREFSEQ_ID\t"
        "DNA_change\tPROTEIN_change\tHOMOZYGOUS\tLOF\t"
        "CLNSIG\tCLNDN\tREVIEW_STATUS\tNEAR_EDIT\tDIST_TO_EDIT\tTIER\n"
    )

    with open(f"{folder}/{vcf_file}") as vcf, \
         open(f"{folder}/{output_file}", "w") as out:

        out.write(header)

        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            pos_str = fields[1]
            quality_str = fields[5]
            info_str = fields[7]
            format_fields = fields[8]
            sample_field = fields[9].strip() if len(fields) > 9 else ""

            try:
                quality = float(quality_str)
            except ValueError:
                quality = 0.0
            qual_scores.append(quality)

            # Genotype
            fmt_keys = format_fields.split(":")
            samp_vals = sample_field.split(":")
            fmt_dict = dict(zip(fmt_keys, samp_vals)) if len(fmt_keys) == len(samp_vals) else {}
            gt = fmt_dict.get("GT", "./.")
            alleles = gt.replace("|", "/").split("/")
            homozygous = len(set(alleles)) == 1

            # Distance to nearest edit site
            near_edit = False
            dist_to_edit = "."
            pos = int(pos_str)
            if chrom in edit_coords:
                dists = [abs(pos - s) for s in edit_coords[chrom]]
                min_dist = min(dists)
                if min_dist <= window_size:
                    near_edit = True
                    dist_to_edit = str(min_dist)

            # Parse INFO
            info = _parse_info(info_str)
            clnsig, clndn, revstat, clnsig_display = _parse_clnsig(info)

            # Determine CLNSIG tier for filtering
            is_pathogenic = clnsig == "Pathogenic_or_LP"
            is_benign = clnsig == "Benign_or_LB"

            # Parse ANN and LOF fields
            ann_field = info.get("ANN", "")
            lof_field = info.get("LOF", "")
            lof_gene = False
            if lof_field:
                try:
                    lof_gene = lof_field.strip("()").split("|")[0]
                    lof_gene = bool(lof_gene)
                except Exception:
                    lof_gene = False

            annotations = ann_field.split(",") if ann_field else []

            for ann in annotations:
                args = ann.split("|")
                if len(args) < 11:
                    continue

                effect_type = args[1]   # e.g. missense_variant
                impact = args[2]        # MODIFIER / LOW / MODERATE / HIGH
                gene = args[3]
                transcript = args[6]
                dna_change = args[9]
                prot_change = args[10]

                impact_idx = intensities.index(impact) if impact in intensities else 0

                # --- Tiering ---
                # Priority: NEAR_EDIT > P/LP > Benign (exclude) > HIGH > MODERATE > rest
                # ClinVar expert curation overrides SnpEff impact predictions for known variants.
                tier = None

                if near_edit:
                    tier = "NEAR_EDIT"
                elif is_pathogenic:
                    tier = "TIER1_CLINVAR_PLP"
                elif is_benign:
                    tier = None  # exclude — ClinVar B/LB supersedes SnpEff impact
                elif impact == "HIGH" or effect_type in _HIGH_EFFECTS:
                    tier = "TIER1_HIGH_IMPACT"
                elif impact == "MODERATE":
                    tier = "TIER2_MODERATE"
                elif impact_idx < min_idx:
                    tier = None  # below minimum intensity → exclude

                # Apply user min_intensity override: if variant is below minimum
                # intensity and not otherwise elevated, skip (unless near edit)
                if tier is None:
                    continue

                # Quality gate (not applied to NEAR_EDIT variants)
                if tier != "NEAR_EDIT" and quality < float(min_quality):
                    continue

                row = [
                    f"{chrom}:{pos_str}",
                    str(quality),
                    effect_type,
                    impact,
                    gene,
                    transcript,
                    dna_change,
                    prot_change,
                    str(homozygous),
                    str(lof_gene),
                    clnsig_display if clnsig_display else ".",
                    clndn if clndn else ".",
                    revstat if revstat else ".",
                    str(near_edit),
                    dist_to_edit,
                    tier,
                ]
                out.write("\t".join(row) + "\n")

    # Quality histogram
    if qual_scores:
        plt.figure()
        plt.hist(qual_scores, bins=[i for i in range(0, 250, 10)])
        plt.title(f"SNP Quality Histogram — {output_file}")
        plt.xlabel("Quality Score")
        plt.ylabel("Frequency")
        plt.savefig(f"{folder}/seqverify_snp_quality.png")
        plt.close()
