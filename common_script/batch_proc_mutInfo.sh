#! /bin/sh
#$ -S /bin/sh
#$ -cwd

INPUT=$1
OUTPUTDIR=$2
LABEL=$3

binary_path=`dirname $0`

echo "perl ${binary_path}/proc_ind.pl ${INPUT} 3 0 > ${OUTPUTDIR}/${LABEL}.ind3.txt"
perl ${binary_path}/proc_ind.pl ${INPUT} 3 0 > ${OUTPUTDIR}/${LABEL}.ind3.txt

echo "perl ${binary_path}/proc_ind.pl ${INPUT} 3 1 > ${OUTPUTDIR}/${LABEL}.ind3_dir.txt"
perl ${binary_path}/proc_ind.pl ${INPUT} 3 1 > ${OUTPUTDIR}/${LABEL}.ind3_dir.txt

echo "perl ${binary_path}/proc_ind.pl ${INPUT} 5 0 > ${OUTPUTDIR}/${LABEL}.ind5.txt"
perl ${binary_path}/proc_ind.pl ${INPUT} 5 0 > ${OUTPUTDIR}/${LABEL}.ind5.txt

echo "perl ${binary_path}/proc_ind.pl ${INPUT} 5 1 > ${OUTPUTDIR}/${LABEL}.ind5_dir.txt"
perl ${binary_path}/proc_ind.pl ${INPUT} 5 1 > ${OUTPUTDIR}/${LABEL}.ind5_dir.txt

echo "perl ${binary_path}/proc_full.pl ${INPUT} 3 0 > ${OUTPUTDIR}/${LABEL}.full3.txt"
perl ${binary_path}/proc_full.pl ${INPUT} 3 0 > ${OUTPUTDIR}/${LABEL}.full3.txt

echo "perl ${binary_path}/proc_full.pl ${INPUT} 3 1 > ${OUTPUTDIR}/${LABEL}.full3_dir.txt"
perl ${binary_path}/proc_full.pl ${INPUT} 3 1 > ${OUTPUTDIR}/${LABEL}.full3_dir.txt

echo "perl ${binary_path}/proc_full.pl ${INPUT} 5 0 > ${OUTPUTDIR}/${LABEL}.full5.txt"
perl ${binary_path}/proc_full.pl ${INPUT} 5 0 > ${OUTPUTDIR}/${LABEL}.full5.txt

echo "perl ${binary_path}/proc_full.pl ${INPUT} 5 1 > ${OUTPUTDIR}/${LABEL}.full5_dir.txt"
perl ${binary_path}/proc_full.pl ${INPUT} 5 1 > ${OUTPUTDIR}/${LABEL}.full5_dir.txt

