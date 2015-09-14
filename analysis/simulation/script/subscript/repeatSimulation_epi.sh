#! /bin/sh
#$ -S /bin/sh
#$ -cwd

# export R_LIBS=/home/yshira/.R
# export R_PATH=/home/yshira/local/bin
export R_PATH=/usr/bin

FEATNUM=$1
SAMPLENUM=$2
ALPHA=$3
GAMMA=$4
OUTPUTFILE=$5

echo "${R_PATH}/R --vanilla --slave --args ${FEATNUM} ${SAMPLENUM} ${ALPHA} ${GAMMA} ${OUTPUTFILE} < repeatSimulation_epi.R"
${R_PATH}/R --vanilla --slave --args ${FEATNUM} ${SAMPLENUM} ${ALPHA} ${GAMMA} ${OUTPUTFILE} < repeatSimulation_epi.R
