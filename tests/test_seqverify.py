"""
Unit + smoke tests for the SeqVerify consolidated-improvements branch.

Runnable two ways:
    pytest tests/                     # if pytest is installed
    python tests/test_seqverify.py    # standalone (no pytest needed)

Covers the pure-Python interpretation modules (CNV filtering, LOH detection,
pathogenicity tiering), module importability, and the CLI/help of the main
script and the plasmid-to-edits helper.
"""
import os
import sys
import subprocess
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


# --------------------------------------------------------------------------- #
def test_imports():
    """All interpretation modules import cleanly."""
    import seqver_pathogenicity          # noqa: F401
    import seqver_liftover               # noqa: F401
    import seqver_loh                    # noqa: F401
    import seqver_cnv_filter             # noqa: F401


# --------------------------------------------------------------------------- #
def test_cnv_filter():
    """filter_cnv_calls assigns PASS / SMALL / REPEAT_MAPQ / NOT_SIG correctly."""
    from seqver_cnv_filter import filter_cnv_calls
    hdr = "type\tchrom\tstart\tend\tsize\tcnv\tp_val\tp_val_2\tp_val_3\tp_val_4\tQ0\tpN\tdG\n"
    rows = [
        # PASS: 500 kb, low Q0, significant
        "deletion\tchr1\t1000000\t1500000\t500000\t0.5\t1e-10\t1e-10\t1e-10\t1e-10\t0.10\t0.0\t1e6",
        # SMALL: 100 kb (< 3*bin_size = 300 kb)
        "duplication\tchr2\t2000000\t2100000\t100000\t1.5\t1e-10\t1e-10\t1e-10\t1e-10\t0.10\t0.0\t2e6",
        # REPEAT_MAPQ: Q0 > 0.5
        "deletion\tchr3\t3000000\t3500000\t500000\t0.5\t1e-10\t1e-10\t1e-10\t1e-10\t0.80\t0.0\t3e6",
        # NOT_SIG: all p-values > 0.05
        "duplication\tchr4\t4000000\t4500000\t500000\t1.5\t0.9\t0.9\t0.9\t0.9\t0.10\t0.0\t4e6",
    ]
    with tempfile.TemporaryDirectory() as d:
        calls = os.path.join(d, "calls.100000.tsv")
        out = os.path.join(d, "calls.100000_filtered.tsv")
        with open(calls, "w") as f:
            f.write(hdr + "\n".join(rows) + "\n")
        filter_cnv_calls(calls, out, bin_size=100000)
        got = {}
        with open(out) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                got[p[1]] = p[-1]          # chrom -> FILTER
        assert got["chr1"] == "PASS", got
        assert got["chr2"] == "SMALL", got
        assert "REPEAT_MAPQ" in got["chr3"], got
        assert "NOT_SIG" in got["chr4"], got
        # PASS-only file has exactly the chr1 call
        pass_path = out.replace(".tsv", "_pass.tsv")
        with open(pass_path) as f:
            body = [l for l in f.read().splitlines()[1:] if l.strip()]
        assert len(body) == 1 and body[0].split("\t")[1] == "chr1", body


# --------------------------------------------------------------------------- #
def _write_loh_vcf(path, chrom, positions, ad_ref, ad_alt):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for p in positions:
            f.write(f"{chrom}\t{p}\t.\tA\tG\t100\t.\t.\tGT:AD\t0/1:{ad_ref},{ad_alt}\n")


def test_loh():
    """detect_loh_around_edits flags balanced vs skewed allele-balance windows."""
    from seqver_loh import detect_loh_around_edits
    with tempfile.TemporaryDirectory() as d:
        center = 5_000_000
        pos = [center + i * 1000 for i in range(8)]     # 8 het SNPs in-window
        # Balanced: AB ~0.5
        vcf_bal = os.path.join(d, "bal.vcf")
        _write_loh_vcf(vcf_bal, "chr1", pos, 10, 10)
        edits = os.path.join(d, "edits.txt")
        with open(edits, "w") as f:
            f.write(f"chr1:{center}-{center}\tACGT\n")
        out_bal = os.path.join(d, "bal_loh.tsv")
        detect_loh_around_edits(vcf_bal, None, out_bal, command_file=edits)
        assert "BALANCED" in open(out_bal).read(), open(out_bal).read()

        # Skewed: AB ~0.9 -> LOH candidate
        vcf_loh = os.path.join(d, "loh.vcf")
        _write_loh_vcf(vcf_loh, "chr1", pos, 2, 18)
        out_loh = os.path.join(d, "loh_loh.tsv")
        detect_loh_around_edits(vcf_loh, None, out_loh, command_file=edits)
        assert "LOH_CANDIDATE" in open(out_loh).read(), open(out_loh).read()


# --------------------------------------------------------------------------- #
def _ann(effect, impact, gene):
    # SnpEff ANN: Allele|Annotation|Impact|Gene|GeneID|Feature|FeatureID|BioType|Rank|HGVS.c|HGVS.p|...
    return f"ANN=G|{effect}|{impact}|{gene}|{gene}|transcript|NM_1.1|protein_coding|1/1|c.1A>G|p.X1Y||||"


def test_pathogenicity_tiers():
    """pathogenicity_logger tiers HIGH / ClinVar-P/LP / MODERATE and drops benign."""
    from seqver_pathogenicity import pathogenicity_logger
    rows = [
        # HIGH-impact stop_gained -> TIER1_HIGH_IMPACT
        ("chr1", 100, "stop", _ann("stop_gained", "HIGH", "GENEH")),
        # ClinVar Pathogenic missense -> TIER1_CLINVAR_PLP
        ("chr2", 200, "path", _ann("missense_variant", "MODERATE", "GENEP") + ";CLNSIG=Pathogenic"),
        # MODERATE missense -> TIER2_MODERATE
        ("chr3", 300, "mod", _ann("missense_variant", "MODERATE", "GENEM")),
        # ClinVar Benign but HIGH ANN -> excluded (benign supersedes)
        ("chr4", 400, "ben", _ann("stop_gained", "HIGH", "GENEB") + ";CLNSIG=Benign"),
    ]
    with tempfile.TemporaryDirectory() as d:
        vcf = "in.vcf"
        with open(os.path.join(d, vcf), "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
            for chrom, pos, _id, info in rows:
                f.write(f"{chrom}\t{pos}\t.\tA\tG\t100\t.\t{info}\tGT:AD\t0/1:5,5\n")
        pathogenicity_logger(d, "out.tsv", vcf, min_quality=3)
        lines = open(os.path.join(d, "out.tsv")).read().splitlines()
        tier = {}
        for l in lines[1:]:
            p = l.split("\t")
            tier.setdefault(p[4], p[-1])       # GENE -> TIER (first row per gene)
        assert tier.get("GENEH") == "TIER1_HIGH_IMPACT", tier
        assert tier.get("GENEP") == "TIER1_CLINVAR_PLP", tier
        assert tier.get("GENEM") == "TIER2_MODERATE", tier
        assert "GENEB" not in tier, f"benign variant should be excluded: {tier}"


# --------------------------------------------------------------------------- #
def test_cli_help():
    """Main script and plasmid helper expose --help and the new toggles."""
    env = dict(os.environ)
    h = subprocess.run([sys.executable, os.path.join(REPO, "seqverify"), "--help"],
                       capture_output=True, text=True, env=env)
    assert h.returncode == 0, h.stderr
    assert "--deduplicate" in h.stdout and "--no-deduplicate" in h.stdout
    assert "--analyze_nonhdr" in h.stdout

    p = subprocess.run([sys.executable, os.path.join(REPO, "scripts", "make_edits_from_plasmid.py"), "--help"],
                       capture_output=True, text=True, env=env)
    assert p.returncode == 0, p.stderr
    assert "--plasmid" in p.stdout and "--keep-rox" in p.stdout


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except Exception as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} tests passed")
    sys.exit(1 if failed else 0)
