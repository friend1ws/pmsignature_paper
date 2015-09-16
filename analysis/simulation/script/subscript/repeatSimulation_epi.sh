#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_PATH=/usr/local/package/r/3.2.1/bin
export R_LIBS=~/.R

FEATNUM=$1
SAMPLENUM=$2
ALPHA=$3
GAMMA=$4
OUTPUTFILE=$5

echo "${R_PATH}/R --vanilla --slave --args ${FEATNUM} ${SAMPLENUM} ${ALPHA} ${GAMMA} ${OUTPUTFILE} < repeatSimulation_epi.R"
${R_PATH}/R --vanilla --slave --args ${FEATNUM} ${SAMPLENUM} ${ALPHA} ${GAMMA} ${OUTPUTFILE} < repeatSimulation_epi.R
