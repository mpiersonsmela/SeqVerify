#!bin/bash
mkdir -p $PREFIX/bin

chmod SeqVerify-prerelease/seqverify +x
cp SeqVerify-prerelease/seqverify $PREFIX/bin
cp SeqVerify-prerelease/seqver_functions.py $PREFIX/bin
cp SeqVerify-prerelease/seqver_transgenes.py $PREFIX/bin
