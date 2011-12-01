#!/bin/bash
# this will require ~ 300 Mb of space

[ -f ../bwa ] || {
	echo "first compile bwa"
	exit;
}

mkdir -p out log/gcov

#Download c.elegans genome fasta:
wget ftp://ftp.ensembl.org/pub/release-64/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.gz

#index fasta: make sure you use the right version of bwa
../bwa index Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.gz

zcat ce_med.fq.gz | head -n50000 | gzip > ce_small.fq.gz


echo you are now ready to run test_diff.sh and test_prof.sh
