#!/bin/bash
# seqver_nonhdr.sh — reconstruct + annotate the non-HDR allele at each targeted
# edit site of a finished SeqVerify sample.
#
# For each edit site it: (1) extracts the read pairs spanning the locus and
# realigns them to a wild-type mini-reference (edit interval +-10 kb); HDR reads
# soft-clip at the cassette junction while non-HDR reads span the cut, (2)
# classifies reads WT vs "scar" and annotates the dominant non-HDR indel for CDS
# impact + distance to the nearest splice site, (3) de-novo assembles the
# non-HDR-only reads with SPAdes and realigns the contigs, and (4) builds a
# standalone IGV report (reads + contig + single-MANE gene model).
#
# All external tools are taken from PATH (samtools, bwa-mem2, spades.py,
# create_report [igv-reports], python (with pysam), passed via --python).
#
# Usage:
#   seqver_nonhdr.sh --bam <insertion.bam> --edits <edits.txt> \
#                    --genome <ref.fa> --gtf <ref.gff3> --out <outdir> \
#                    [--threads N]
#
# Edits file format (one site per line, tab-separated):
#   CHR:START-END<TAB>INSERTED_SEQUENCE
set -u

BAM="" ; EDITS="" ; REF="" ; GFF="" ; OUT="" ; T=8 ; PY="python"
while [ $# -gt 0 ]; do
  case "$1" in
    --bam)     BAM="$2"; shift 2 ;;
    --edits)   EDITS="$2"; shift 2 ;;
    --genome)  REF="$2"; shift 2 ;;
    --gtf)     GFF="$2"; shift 2 ;;
    --out)     OUT="$2"; shift 2 ;;
    --threads) T="$2"; shift 2 ;;
    --python)  PY="$2"; shift 2 ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

for v in BAM EDITS REF GFF OUT; do
  eval "val=\${$v}"
  [ -n "$val" ] || { echo "Missing required --${v,,}"; exit 2; }
done
[ -f "$EDITS" ] || { echo "NO edits file: $EDITS"; exit 1; }
[ -f "$BAM" ]   || { echo "NO bam file: $BAM"; exit 1; }
[ -f "$REF" ]   || { echo "NO genome fasta: $REF"; exit 1; }
for tool in samtools bwa-mem2 spades.py create_report "$PY"; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Required tool not on PATH: $tool"; exit 3; }
done

mkdir -p "$OUT"
[ -f "$REF.fai" ] || samtools faidx "$REF"
[ -f "$BAM.bai" ] || [ -f "$BAM.csi" ] || samtools index -@"$T" "$BAM"
SUM=$OUT/SUMMARY.txt
echo "non-HDR allele analysis  ($(date))"  > "$SUM"
echo "extract locus read-pairs -> realign to WT (+-10kb); HDR reads clip at the cassette junction, non-HDR reads span the cut." >> "$SUM"; echo >> "$SUM"

while IFS=$'\t' read -r loc seq; do
  [ -z "$loc" ] && continue
  chrom=${loc%%:*}; rng=${loc##*:}; START=${rng%-*}; END=${rng#*-}
  sid="${chrom}_${START}"; win_start=$((START-10000)); win_end=$((END+10000)); ext_end=$((END+12000))
  [ "$win_start" -lt 1 ] && win_start=1
  ref=$OUT/${sid}_wt.fa
  samtools faidx "$REF" ${chrom}:${win_start}-${win_end} | sed "1s/.*/>${sid}_win/" > "$ref"; samtools faidx "$ref"
  samtools view "$BAM" ${chrom}:${win_start}-${ext_end} | cut -f1 | sort -u > $OUT/${sid}_names.txt
  samtools view -@"$T" -N $OUT/${sid}_names.txt -b "$BAM" -o $OUT/${sid}_pairs.bam
  samtools collate -@"$T" -O -u $OUT/${sid}_pairs.bam | samtools fastq -@"$T" -1 $OUT/${sid}_R1.fastq -2 $OUT/${sid}_R2.fastq -0 /dev/null -s /dev/null -n 2>/dev/null
  bwa-mem2 index "$ref" >/dev/null 2>&1
  bwa-mem2 mem -t "$T" "$ref" $OUT/${sid}_R1.fastq $OUT/${sid}_R2.fastq 2>/dev/null | samtools sort -@"$T" -o $OUT/${sid}_wtaln.bam - 2>/dev/null
  samtools index $OUT/${sid}_wtaln.bam
  awk -F'\t' -v c=$chrom -v s=$win_start -v e=$win_end '$1==c && $4<=e && $5>=s && ($3=="exon"||$3=="CDS"||$3=="gene"||$3=="transcript"||$3=="mRNA")' "$GFF" > $OUT/${sid}_gff.txt
  ctg="${sid}_win"; istart=$((START-win_start+1)); iend=$((END-win_start+1))
  "$PY" - "$OUT/${sid}_wtaln.bam" "$ctg" "$istart" "$iend" "$win_start" "$chrom" "$sid" "$SUM" "$START" "$END" "$OUT/${sid}_gff.txt" <<'PY'
import pysam,sys,re
bam,ctg,istart,iend,off,chrom,sid,sumf,gS,gE,gff=sys.argv[1:12]
istart=int(istart);iend=int(iend);off=int(off);gS=int(gS);gE=int(gE)
b=pysam.AlignmentFile(bam); lo=istart-25; hi=iend+25
def reads():
    for r in b.fetch(ctg, max(0,lo-60), hi+60):
        if not (r.is_unmapped or r.is_secondary or r.is_supplementary): yield r
# PASS1: catalog indels in interval+-25
indels={}
for r in reads():
    pos=r.reference_start; qpos=0
    for op,ln in (r.cigartuples or []):
        if op in (0,7,8): pos+=ln; qpos+=ln
        elif op==1:
            if lo<=pos<=hi:
                k=f"+{ln}bp:{r.query_sequence[qpos:qpos+ln]}"; indels.setdefault(k,[0,off+pos-1,ln,'ins',pos]); indels[k][0]+=1
            qpos+=ln
        elif op==2:
            if lo<=pos<=hi:
                k=f"-{ln}bp"; indels.setdefault(k,[0,off+pos,ln,'del',pos]); indels[k][0]+=1
            pos+=ln
        elif op==4: qpos+=ln
dom=sorted(indels.items(),key=lambda x:-x[1][0])
P = dom[0][1][4] if dom else (istart+iend)//2   # scar refpos (0-based window)
# PASS2: classify reads spanning the cut P
wt=scar=clip=part=0
for r in reads():
    cig=r.cigartuples or []; clipped=False
    if cig:
        if cig[0][0]==4 and (P-25)<=r.reference_start<=(P+25): clipped=True
        if cig[-1][0]==4 and (P-25)<=r.reference_end<=(P+25): clipped=True
    if not (r.reference_start<P-5 and r.reference_end>P+5):
        clip+=(1 if clipped else 0); part+=(0 if clipped else 1); continue
    pos=r.reference_start; has=False
    for op,ln in cig:
        if op in (0,7,8): pos+=ln
        elif op==1:
            if P-10<=pos<=P+10: has=True
        elif op==2:
            if P-10<=pos<=P+10: has=True
            pos+=ln
    if has: scar+=1
    else: wt+=1
def parse_gff():
    genes=[];feats=[]
    for l in open(gff):
        f=l.rstrip("\n").split("\t")
        if len(f)<9: continue
        t,s,e,st,at=f[2],int(f[3]),int(f[4]),f[6],f[8]
        if t=="gene":
            gn=re.search(r'gene_name=([^;]+)',at) or re.search(r'gene=([^;]+)',at); genes.append((s,e,st,gn.group(1) if gn else "?"))
        elif t in ("exon","CDS"):
            par=re.search(r'Parent=([^;]+)',at); gg=re.search(r'gene=([^;]+)',at)
            feats.append((t,s,e,st,par.group(1) if par else "?",gg.group(1) if gg else "?", 'MANE Select' in at))
    return genes,feats
def annotate(gpos,isize,itype):
    genes,feats=parse_gff()
    gc=[g for g in genes if g[0]<=gE and g[1]>=gS]; gene=gc[0][3] if gc else "?"
    txs={}
    for t,s,e,st,par,gg,mane in feats:
        if gg!=gene: continue
        d=txs.setdefault(par,{'exons':[],'cds':[],'mane':False,'strand':st})
        (d['exons'] if t=="exon" else d['cds']).append((s,e))
        if mane: d['mane']=True
    if not txs: return f"gene={gene}; no transcript in window"
    txid,txd=sorted(txs.items(),key=lambda kv:(not kv[1]['mane'],-len(kv[1]['cds'])))[0]
    ex=sorted(txd['exons']); cds=sorted(txd['cds']); strand=txd['strand']
    in_cds=any(s<=gpos<=e for s,e in cds); in_exon=any(s<=gpos<=e for s,e in ex)
    bnds=[]
    for i,(s,e) in enumerate(ex):
        if i>0: bnds.append(('acceptor' if strand=='+' else 'donor',s))
        if i<len(ex)-1: bnds.append(('donor' if strand=='+' else 'acceptor',e))
    if bnds:
        nb=min(bnds,key=lambda x:abs(x[1]-gpos)); dist=abs(nb[1]-gpos)
        sstxt=(f"{dist} bp from nearest splice site ({nb[0]} at {chrom}:{nb[1]})" if dist<=100 else "no splice site within 100 bp")
    else: sstxt="single-exon transcript (no splice sites)"
    cds_txt=("IN CDS — "+("FRAMESHIFT" if isize%3 else f"in-frame {itype} {isize}bp")) if in_cds else ("NOT in CDS ("+("exonic/UTR" if in_exon else "intronic/intergenic")+")")
    return f"gene={gene} tx={txid} ({strand}); {cds_txt}; {sstxt}"
with open(sumf,"a") as o:
    o.write(f"SITE {sid}  edit interval {chrom}:{gS}-{gE}\n")
    if dom:
        d=dom[0]; cnt,gpos,isize,itype=d[1][0],d[1][1],d[1][2],d[1][3]
        o.write(f"  non-HDR scar: {d[0]} at {chrom}:{gpos}  ({cnt} reads carry it)\n")
        o.write(f"  annotation: {annotate(gpos,isize,itype)}\n")
        if len(dom)>1: o.write(f"  other indels nearby: {{{', '.join(k+':'+str(v[0]) for k,v in dom[1:] if v[0]>1)}}}\n")
    else: o.write("  non-HDR scar: NONE detected\n")
    o.write(f"  spanning cut: WT={wt}, scar={scar}; HDR-clipped={clip}; partial={part}\n")
    if clip>=5 and scar>=3 and wt==0: o.write("  => HET: HDR / non-HDR indel (no WT allele)\n")
    elif clip>=5 and wt>=3 and scar==0: o.write("  => HET: HDR / WILD-TYPE non-HDR allele\n")
    elif clip<3 and scar>=3: o.write("  => HOMOZYGOUS indel or HDR/indel (few junction reads)\n")
    else: o.write("  => ambiguous — inspect\n")
PY
  samtools view $OUT/${sid}_wtaln.bam | awk '$6 ~ /([3-9][0-9]|[0-9]{3,})S/ {print $1}' > $OUT/${sid}_excl.txt
  samtools view -f4 $OUT/${sid}_wtaln.bam 2>/dev/null | cut -f1 >> $OUT/${sid}_excl.txt; sort -u $OUT/${sid}_excl.txt -o $OUT/${sid}_excl.txt
  "$PY" - "$OUT" "$sid" <<'PY'
import sys; OUT,sid=sys.argv[1:3]; excl=set(open(f"{OUT}/{sid}_excl.txt").read().split())
for rd in ["R1","R2"]:
    with open(f"{OUT}/{sid}_{rd}.fastq") as f,open(f"{OUT}/{sid}_nonHDR_{rd}.fastq","w") as o:
        while True:
            h=f.readline()
            if not h: break
            s=f.readline();p=f.readline();q=f.readline()
            if h[1:].split()[0] not in excl: o.write(h+s+p+q)
PY
  rm -rf $OUT/${sid}_spades
  spades.py --careful -1 $OUT/${sid}_nonHDR_R1.fastq -2 $OUT/${sid}_nonHDR_R2.fastq -o $OUT/${sid}_spades -t "$T" >$OUT/${sid}_spades.log 2>&1
  if [ -s $OUT/${sid}_spades/contigs.fasta ]; then
    bwa-mem2 mem -t "$T" "$ref" $OUT/${sid}_spades/contigs.fasta 2>/dev/null | samtools sort -@"$T" -o $OUT/${sid}_ctg_vs_wt.bam - 2>/dev/null
    samtools index $OUT/${sid}_ctg_vs_wt.bam; mid=$(((istart+iend)/2))
    cinfo=$(samtools view $OUT/${sid}_ctg_vs_wt.bam ${ctg}:$((mid-3))-$((mid+3)) | awk 'NR==1{print $1" CIGAR "substr($6,1,90)}')
    echo "  assembly: ${cinfo:-no contig spans cut}" >> "$SUM"
  else echo "  assembly: no contigs" >> "$SUM"; fi
  # ---- IGV viewer for the non-HDR allele (reads realigned to WT + assembled contig + gene model) ----
  printf "%s\t%d\t%d\n" "$ctg" $((istart-1)) "$iend" > $OUT/${sid}_view.bed
  WL=$((win_end-win_start+1))
  # clean single MANE-transcript gene model (mRNA+exon+CDS) so igv.js renders the CDS boundary
  "$PY" - "$OUT/${sid}_gff.txt" "$ctg" $((win_start-1)) "$WL" "$OUT/${sid}_track.gff3" <<'PY'
import sys,re
from collections import Counter
gff,ctg,off,wl,out=sys.argv[1:6]; off=int(off); wl=int(wl)
L=[l.rstrip("\n").split("\t") for l in open(gff) if not l.startswith("#") and len(l.rstrip("\n").split("\t"))>=9]
txspan={}; gene_of={}; mane=None
for f in L:
    if f[2] in ("transcript","mRNA"):
        m=re.search(r'ID=([^;]+)',f[8])
        if not m: continue
        tid=m.group(1); txspan[tid]=(int(f[3]),int(f[4]),f[6])
        g=re.search(r'gene=([^;]+)',f[8]); gene_of[tid]=g.group(1) if g else tid
        if 'MANE Select' in f[8] and mane is None: mane=tid
if mane is None:
    c=Counter()
    for f in L:
        if f[2]=="CDS":
            m=re.search(r'Parent=([^;]+)',f[8])
            if m: c[m.group(1)]+=1
    mane=c.most_common(1)[0][0] if c else (next(iter(txspan)) if txspan else None)
def rm(s,e):
    s-=off; e-=off; return max(1,s),min(wl,e)
with open(out,"w") as o:
    o.write("##gff-version 3\n")
    if mane and mane in txspan:
        s,e,st=txspan[mane]; rs,re_=rm(s,e)
        o.write(f"{ctg}\tLiftoff\tmRNA\t{rs}\t{re_}\t.\t{st}\t.\tID={mane};Name={gene_of.get(mane,mane)}\n")
        feats=[]
        for f in L:
            if f[2] in ("exon","CDS"):
                m=re.search(r'Parent=([^;]+)',f[8])
                if m and m.group(1)==mane:
                    a,b=rm(int(f[3]),int(f[4])); feats.append((a,b,f[2],f[7]))
        for a,b,t,ph in sorted(feats):
            o.write(f"{ctg}\tLiftoff\t{t}\t{a}\t{b}\t.\t{st}\t{ph if t=='CDS' else '.'}\tParent={mane}\n")
PY
  tracks="$OUT/${sid}_wtaln.bam"
  [ -f "$OUT/${sid}_ctg_vs_wt.bam" ] && tracks="$tracks $OUT/${sid}_ctg_vs_wt.bam"
  [ -s "$OUT/${sid}_track.gff3" ] && tracks="$tracks $OUT/${sid}_track.gff3"
  create_report $OUT/${sid}_view.bed --fasta "$ref" --standalone --flanking 400 \
    --sequence 1 --begin 2 --end 3 --tracks $tracks \
    --output $OUT/${sid}_nonHDR_igv.html >$OUT/${sid}_igv.log 2>&1 \
    && echo "  IGV viewer: $OUT/${sid}_nonHDR_igv.html" >> "$SUM" \
    || echo "  IGV viewer: FAILED (see ${sid}_igv.log)" >> "$SUM"
  echo >> "$SUM"
  rm -f $OUT/${sid}_pairs.bam $OUT/${sid}_wt.fa.0123 $OUT/${sid}_wt.fa.bwt.2bit.64 $OUT/${sid}_wt.fa.amb $OUT/${sid}_wt.fa.ann $OUT/${sid}_wt.fa.pac
done < "$EDITS"
echo "DONE non-HDR analysis"; cat "$SUM"
