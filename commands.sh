#!/bin/bash

# Michael E. Sparks (michael.sparks2@usda.gov), 20 April 2023

if [ -n "$1" ]; then
  # $handle will correspond to ``LdIV1.polyprotein" for sample data
  handle=$1
else
  echo "Please supply filename handle as first argument (additional arguments are ignored)"
  exit 1
fi

#$ muscle -version
#MUSCLE v3.8.1551 by Robert C. Edgar
cat ${handle}_PEP.fpa | muscle -clwstrict > ${handle}_PEP.msa

#$ seqret --version
#EMBOSS:6.6.0.0
seqret -sequence ${handle}_PEP.msa -outseq ${handle}_PEP.faln -osformat fasta

cat ${handle}_PEP.faln | ./util/merge_fasta.pl > ${handle}_PEP.faln.4tranalign

grep '^>' ${handle}_PEP.faln.4tranalign | sed 's/^>//' | \
  ./util/extract_fasta_seqs.pl <(cat ${handle}_CDS.fna | ./util/merge_fasta.pl) \
  > ${handle}_CDS.fna.4tranalign

# EMBOSS' tranalign is a clone of Bill Pearson's MRTRANS program.
# I've not tested whether MRTRANS' output is interchangeable with tranalign's here
#$ tranalign --version
#EMBOSS:6.6.0.0
tranalign ${handle}_CDS.fna.4tranalign ${handle}_PEP.faln.4tranalign ${handle}_CDS.faln

seqret -sequence ${handle}_CDS.faln -outseq ${handle}_CDS.nex -osformat nexus

cat ${handle}_PEP.faln | ./util/merge_fasta.pl > ${handle}_PEP.faln.mer
cat ${handle}_CDS.faln | ./util/merge_fasta.pl > ${handle}_CDS.faln.mer

# plain text listing of mutations, w/ synonymous muts flagged by '=', others by '^'
./util/merge_alns.pl ${handle}_PEP.faln.mer ${handle}_CDS.faln.mer | unix2dos > ${handle}_PEPandCDS.txt
# In principle, the following could change if muscle's msa turns out differently from my results.
#$ md5sum LdIV1.polyprotein_PEPandCDS.txt
#9d7421a944b8b311bdf2a01b71325948  LdIV1.polyprotein_PEPandCDS.txt

# tabulates mutations for slurping into Excel, if desired
./util/tab_muts.pl ${handle}_PEP.faln.mer ${handle}_CDS.faln.mer | unix2dos > ${handle}_mut_table.4excel.txt
# ditto on variability w/ muscle output
#$ md5sum LdIV1.polyprotein_mut_table.4excel.txt
#ed961a999617ad4a65e1780831ca62db  LdIV1.polyprotein_mut_table.4excel.txt

echo "Done. Please inspect ${handle}_PEPandCDS.txt and ${handle}_mut_table.4excel.txt"

exit 0
