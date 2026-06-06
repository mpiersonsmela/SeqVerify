"""
seqver_liftover.py
Coordinate adjustment (edited-CHM13 → unedited-CHM13) and CrossMap VCF liftover.

Two-step liftover:
  1. adjust_vcf_for_edits(): rewrite VCF positions from edited genome back to
     unedited-CHM13 coordinates using the known edit list.
  2. liftover_vcf(): call CrossMap to lift unedited-CHM13 VCF → GRCh38.

Wild-type path (no edits): adjust_vcf_for_edits() is a no-op copy.
"""

import os
import sys


def adjust_vcf_for_edits(vcf_in, edits, vcf_out, edit_region_out=None):
    """
    Adjust VCF positions from edited-CHM13 space back to unedited-CHM13 space.

    edits: list of editArgument objects as returned by
           commandHandler(..., return_commands_only=True).
           Each has .chr (str) and .coords [start, end] (0-based Python coords
           in cumulative-edited-genome space, as set by commandHandler).
           .sequence holds the inserted/replacement sequence.

    Variants that fall strictly inside an inserted sequence have no unedited
    equivalent; they are written to edit_region_out (if given) and omitted
    from vcf_out.
    """
    # Build per-chromosome edit list: each entry is
    # (edited_start, edited_end, delta) where delta = len(seq)-(orig_end-orig_start)
    # For the reverse transform we need the positions *in edited-genome space*
    # and the shift each edit introduced.
    #
    # commandHandler sets coords so that each edit's [start,end] is valid in the
    # genome AFTER all previous edits on the same chromosome are applied.
    # We process edits in position order (they are already sorted by commandHandler).

    edits_by_chr = {}
    for cmd in (edits or []):
        chrom = cmd.chr
        seq_len = len(cmd.sequence)
        orig_span = cmd.coords[1] - cmd.coords[0]  # original length replaced
        delta = seq_len - orig_span                 # positive = insertion, negative = deletion
        entry = (cmd.coords[0], cmd.coords[0] + seq_len, delta, orig_span)
        # (edited_start, edited_end_of_new_seq, delta, original_span)
        edits_by_chr.setdefault(chrom, []).append(entry)

    skipped = 0
    kept = 0

    with open(vcf_in) as fin, open(vcf_out, "w") as fout:
        edit_fh = open(edit_region_out, "w") if edit_region_out else None
        try:
            for line in fin:
                if line.startswith("#"):
                    fout.write(line)
                    if edit_fh:
                        edit_fh.write(line)
                    continue

                fields = line.rstrip("\n").split("\t")
                chrom = fields[0]
                # VCF positions are 1-based
                pos = int(fields[1])
                pos0 = pos - 1  # convert to 0-based for arithmetic

                chr_edits = edits_by_chr.get(chrom, [])
                in_edit_region = False
                shift = 0

                for edited_start, edited_end, delta, orig_span in chr_edits:
                    if pos0 < edited_start + shift:
                        # This position is before this edit; no more adjustments.
                        # (edits are ordered, so subsequent ones are further right)
                        break
                    if pos0 < edited_end + shift:
                        # Position falls inside the inserted/replaced sequence —
                        # it has no equivalent in the unedited genome.
                        in_edit_region = True
                        break
                    # Position is after this edit's inserted block.
                    shift -= delta  # undo the edit's coordinate shift

                if in_edit_region:
                    if edit_fh:
                        fields[7] = fields[7] + ";EDIT_REGION" if fields[7] != "." else "EDIT_REGION"
                        edit_fh.write("\t".join(fields) + "\n")
                    skipped += 1
                    continue

                adjusted_pos = pos0 + shift + 1  # back to 1-based
                fields[1] = str(adjusted_pos)
                fout.write("\t".join(fields) + "\n")
                kept += 1
        finally:
            if edit_fh:
                edit_fh.close()

    print(f"adjust_vcf_for_edits: {kept} variants adjusted, {skipped} in edit regions")
    return kept, skipped


def liftover_vcf(vcf_in, chain_file, ref_fasta, vcf_out, unmapped_out=None):
    """
    Liftover a VCF from CHM13 coordinates to GRCh38 using CrossMap CLI.

    CrossMap 0.7.x is a CLI tool; this function calls it via subprocess.
    Command: CrossMap vcf <chain> <in.vcf> <ref.fa> <out.vcf>

    chain_file: path to chm13v2.0ToHg38 (hs1ToHg38.over.chain.gz)
    ref_fasta: GRCh38 reference FASTA (required by CrossMap for allele checking)
    unmapped_out: ignored (CrossMap writes <out_vcf>.unmap automatically)
    """
    import subprocess
    import shutil

    crossmap_bin = shutil.which("CrossMap")
    if crossmap_bin is None:
        # Look in common conda env location
        env_bin = os.path.expanduser(
            "~/miniconda3/envs/seqverify/bin/CrossMap")
        if os.path.isfile(env_bin):
            crossmap_bin = env_bin
        else:
            sys.exit("CrossMap not found. Run: pip install CrossMap")

    if not os.path.exists(chain_file):
        sys.exit(f"Chain file not found: {chain_file}\n"
                 "Download with --download_defaults or set chain_file= in config.")

    # CrossMap requires a ref FASTA; use /dev/null as a dummy if not provided
    # (allele validation will be skipped but coordinates will still be lifted)
    fa_arg = ref_fasta if ref_fasta and os.path.isfile(ref_fasta) else "/dev/null"

    cmd = [crossmap_bin, "vcf", chain_file, vcf_in, fa_arg, vcf_out]
    print(f"CrossMap liftover: {vcf_in} → {vcf_out}")
    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())

    if os.path.exists(vcf_out):
        n = sum(1 for l in open(vcf_out) if not l.startswith("#"))
        print(f"CrossMap: {n} variants lifted to GRCh38")
    else:
        print("Warning: CrossMap produced no output file")
