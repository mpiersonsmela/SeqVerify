#!/usr/bin/env bash
# setup.sh — fetch the extra reference resource that the variant-interpretation
# stage needs but that is not a conda package: the CHM13 -> GRCh38 (hg38)
# liftover chain file. All software dependencies (samtools, bwa-mem2, spades,
# crossmap, snpeff/snpsift on openjdk 21, cnvpytor, igv-reports, pysam, ...) are
# installed by the conda environment (see seqverify-env.yml).
#
# Usage:
#   bash setup.sh [DEST_DIR]        # DEST_DIR defaults to ./seqverify_defaults
set -euo pipefail

DEST="${1:-seqverify_defaults}"
mkdir -p "$DEST"

CHAIN="$DEST/chm13v2.0ToHg38.over.chain.gz"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz"

if [ -s "$CHAIN" ]; then
  echo "Liftover chain already present: $CHAIN"
else
  echo "Downloading CHM13->hg38 liftover chain to $CHAIN ..."
  curl -L -o "$CHAIN" "$CHAIN_URL"
  echo "Done: $CHAIN"
fi

echo ""
echo "Reminder: create/activate the conda environment first, e.g."
echo "  conda env create -f seqverify-env.yml   # or: conda install -c bioconda seqverify"
echo "  conda activate seqverify"
echo ""
echo "The SnpEff GRCh38 MANE database (variant_db, e.g. GRCh38.mane.1.2.refseq) is"
echo "downloaded automatically by SnpEff on first use, or with:"
echo "  java -jar \$(dirname \$(command -v snpEff))/../share/snpeff-*/snpEff.jar download GRCh38.mane.1.2.refseq"
