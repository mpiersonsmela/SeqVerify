#!/usr/bin/env python3
"""
make_edits_from_plasmid.py

Generate a SeqVerify --targeted edits file from an HDR donor plasmid (GenBank).

Method
------
1. Read the plasmid sequence from a GenBank (.gb) file.
2. BLAST the whole plasmid against a reference genome (e.g. T2T-CHM13v2.0).
   The two homology arms show up as the two strongest HSPs at the target
   locus (long, ~100% identity). The cargo between them is the knock-in.
3. Build the edit in SeqVerify command format:  CHR:START-END<TAB>SEQUENCE
      START = genomic coordinate of the last base of the upstream (5') arm
              (this base is kept)
      END   = (genomic start of the downstream (3') arm) - 1
              (bases START+1..END are deleted)
      SEQUENCE = plasmid cargo between the arms, written on the genome (+) strand
   See the SeqVerify README: "START-END is an interval of coordinates to be
   deleted, not including START. SEQUENCE is the sequence to be inserted from
   START onwards."
4. (Default) Collapse the Rox selection cassette: if two Rox sites are present
   as a direct repeat inside the cargo, delete the sequence between them plus
   one Rox site, leaving a single Rox scar -- matching cell lines in which the
   Rox-flanked PGK-PuroTK cassette was excised by Dre before sequencing.
   Disable with --keep-rox.

Usage
-----
  python make_edits_from_plasmid.py \
      --plasmid edits/synapsis_edits/phdr2-six6os1-mscarlet.gb \
      --ref     seqverify_defaults/chm13v2.0.fa \
      --out     edits/SIX6OS1-mScarlet.txt

  # keep the selection cassette (do not excise Rox):
  python ... --keep-rox
"""
import argparse, os, re, subprocess, sys, tempfile

# Dre recombinase target (rox) site. Searched in both orientations.
DEFAULT_ROX = "TAACTTTAAATAATGCCAATTATTTAAAGTTA"

COMP = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
def revcomp(s):
    return "".join(COMP[b] for b in reversed(s))

def read_genbank_seq(path):
    seq, in_origin = [], False
    with open(path) as f:
        for line in f:
            if line.startswith("ORIGIN"):
                in_origin = True; continue
            if in_origin:
                if line.startswith("//"): break
                seq.append(re.sub(r'[^acgtnACGTN]', '', line))
    s = "".join(seq).upper()
    if not s:
        sys.exit(f"ERROR: no ORIGIN sequence found in {path}")
    return s

def ensure_blastdb(ref):
    """Return a blast db prefix for ref, building one in a temp dir if needed."""
    for ext in (".nin", ".nal"):
        if os.path.exists(ref + ext):
            return ref
    cache = os.path.join(tempfile.gettempdir(),
                         "blastdb_" + re.sub(r'\W', '_', os.path.abspath(ref)))
    if not (os.path.exists(cache + ".nin") or os.path.exists(cache + ".nal")):
        os.makedirs(os.path.dirname(cache) or ".", exist_ok=True)
        sys.stderr.write(f"[info] building BLAST db for {ref} ...\n")
        subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "nucl",
                        "-out", cache], check=True,
                       stdout=subprocess.DEVNULL)
    return cache

def blast_arms(plasmid_fa, db, min_arm):
    """Return the two homology-arm HSPs as dicts, plus chrom/strand.

    Picks, among HSPs >= min_arm, the two longest that share one (chrom, strand)
    -- these are the two homology arms flanking the cargo.
    """
    fmt = "6 qseqid sseqid pident length qstart qend sstart send sstrand bitscore"
    out = subprocess.run(
        ["blastn", "-query", plasmid_fa, "-db", db, "-outfmt", fmt,
         "-evalue", "1e-20", "-max_target_seqs", "50"],
        check=True, capture_output=True, text=True).stdout
    hsps = []
    for ln in out.splitlines():
        f = ln.split("\t")
        hsps.append(dict(chrom=f[1], pident=float(f[2]), length=int(f[3]),
                         qs=int(f[4]), qe=int(f[5]), ss=int(f[6]), se=int(f[7]),
                         strand=f[8], bits=float(f[9])))
    hsps = [h for h in hsps if h["length"] >= min_arm]
    if len(hsps) < 2:
        sys.exit("ERROR: fewer than two arm-length HSPs found; "
                 "lower --min-arm or check the reference.")
    # group by (chrom, strand); choose the group whose two longest HSPs are largest
    groups = {}
    for h in hsps:
        groups.setdefault((h["chrom"], h["strand"]), []).append(h)
    best = None
    for key, g in groups.items():
        g.sort(key=lambda x: x["length"], reverse=True)
        if len(g) < 2:
            continue
        score = g[0]["length"] + g[1]["length"]
        if best is None or score > best[0]:
            best = (score, key, g[:2])
    if best is None:
        sys.exit("ERROR: could not find two arms sharing a chromosome/strand.")
    _, (chrom, strand), (a, b) = best
    return chrom, strand, a, b

def find_rox_pair(plasmid, rox):
    """Return (rox1_start, rox1_end, rox2_start, rox2_end) 1-based, or None.

    Requires exactly two direct-repeat occurrences (same orientation).
    """
    for motif in (rox, revcomp(rox)):
        hits = [m.start() + 1 for m in re.finditer(motif, plasmid)]
        if len(hits) == 2:
            s1, s2 = sorted(hits)
            L = len(motif)
            return (s1, s1 + L - 1, s2, s2 + L - 1)
    return None

def build_edit(plasmid, chrom, strand, a, b, remove_rox, rox):
    # genome (+) strand interval of each arm
    def ginterval(h): return (min(h["ss"], h["se"]), max(h["ss"], h["se"]))
    ia, ib = ginterval(a), ginterval(b)
    armL, armU = (ia, ib) if ia[0] < ib[0] else (ib, ia)
    START = armL[1]          # last base of upstream arm (kept)
    END   = armU[0] - 1      # last deleted base (upstream of downstream arm)

    # cargo = plasmid between the two HSPs (1-based inclusive)
    first, second = (a, b) if a["qs"] < b["qs"] else (b, a)
    c1, c2 = first["qe"] + 1, second["qs"] - 1     # cargo plasmid coords

    # indices (1-based) of cargo to keep
    keep = set(range(c1, c2 + 1))
    rox_note = "kept"
    if remove_rox:
        pair = find_rox_pair(plasmid, rox)
        if pair:
            r1s, r1e, r2s, r2e = pair
            # excise everything after the first rox through the second rox,
            # leaving a single rox scar
            for i in range(r1e + 1, r2e + 1):
                keep.discard(i)
            rox_note = f"excised {r1e+1}-{r2e} (single rox scar retained)"
        else:
            rox_note = "no direct-repeat rox pair found; nothing excised"

    cargo = "".join(plasmid[i - 1] for i in range(c1, c2 + 1) if i in keep)
    insert = revcomp(cargo) if strand == "minus" else cargo
    return START, END, insert, rox_note

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plasmid", required=True, help="GenBank (.gb) donor plasmid")
    ap.add_argument("--ref", required=True, help="reference genome FASTA (e.g. chm13v2.0.fa)")
    ap.add_argument("--out", required=True, help="output edits .txt file")
    ap.add_argument("--name", default=None, help="label for logging (default: output basename)")
    ap.add_argument("--min-arm", type=int, default=150, help="min HSP length for an arm (bp)")
    ap.add_argument("--rox", default=DEFAULT_ROX, help="rox site sequence")
    ap.add_argument("--keep-rox", action="store_true",
                    help="do NOT excise the Rox-flanked selection cassette")
    args = ap.parse_args()
    name = args.name or os.path.splitext(os.path.basename(args.out))[0]

    plasmid = read_genbank_seq(args.plasmid)
    db = ensure_blastdb(args.ref)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        tf.write(f">{name}\n{plasmid}\n"); qfa = tf.name
    try:
        chrom, strand, a, b = blast_arms(qfa, db, args.min_arm)
    finally:
        os.unlink(qfa)

    START, END, insert, rox_note = build_edit(
        plasmid, chrom, strand, a, b, remove_rox=not args.keep_rox, rox=args.rox)

    with open(args.out, "w") as o:
        o.write(f"{chrom}:{START}-{END}\t{insert}\n")

    sys.stderr.write(
        f"[{name}] {chrom}:{START}-{END}  strand={strand}  "
        f"deletion={END-START}bp  insert={len(insert)}bp  rox={rox_note}\n"
        f"          arms: {a['length']}bp ({a['pident']}%) + {b['length']}bp ({b['pident']}%)\n"
        f"          wrote {args.out}\n")

if __name__ == "__main__":
    main()
