#!/bin/bash
seqkit grep -f ${1}_goodhit.txt ${1}_50mer.fasta > ${1}_indexes.fasta
makeblastdb -in ${1}.fasta -dbtype nucl -parse_seqids
megablast -d ${1}.fasta -i ${1}_indexes.fasta -m 8 | sort -k 9 -n > ${1}.tab
mkdir OUTPUT
cp ${1}.tab OUTPUT/

