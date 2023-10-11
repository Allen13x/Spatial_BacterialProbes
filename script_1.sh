#!/bin/bash

sed -z 's/,/\n/g' ${1}50mer |sed 's/./&-/25' | sort -u > ${1}50mer_sep
grep -i "\-A" ${1}50mer_sep  | grep -i -v -P "AAAA|TTTT|CCCC|GGGG" > ${1}50mer_sep_hom3
paste ${1}50mer_sep_hom3 <(cat ${1}50mer_sep_hom3 | tr -d '-' | tr -d -c 'cgCG\n' | awk '{print length/50*100}') | awk -F '\t' ' $2 > 50 && $2 < 70 ' > ${1}50mer_sep_hom3_gc50-70
i=1
while read p; do echo ">k_${i}"; echo $(echo ${p} | tr -d '-' | cut -f1 -d ' '); ((i = i + 1)); done < ${1}50mer_sep_hom3_gc50-70 > ${1}_50mer.fasta

