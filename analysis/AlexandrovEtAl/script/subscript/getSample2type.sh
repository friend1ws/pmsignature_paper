#! /bin/sh

# command for obtaining the relationship between sample and cancer types
while read type; do for sample in `zcat ../../result/MPFormat/${type}.mp.txt.gz | cut -f 1 | sort -u`; do echo -e "${sample}\t${type}"; done; done < ../../data/cancer_types.txt > ../../../figure/data/AlexandrovEtAl_sample2type.txt

