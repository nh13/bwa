#!/bin/bash

FASTA=e_coli_dh10b.fasta;
if [ ! -f ${FASTA} ]; then
	gunzip -c ${FASTA}.gz > ${FASTA};
fi
if [ ! -f ${FASTA}.sa ]; then
	../bwa index ${FASTA};
fi
../bwa speed -f ${FASTA};
