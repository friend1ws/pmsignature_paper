#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

MUTNUM=$1
SAMPLENUM=$2
ALPHA=$3
GAMMA=$4
OUTPUTFILE=$5

echo "${R_PATH}/R --vanilla --slave --args ${MUTNUM} ${SAMPLENUM} ${ALPHA} ${GAMMA} ${OUTPUTFILE} < repeatSimulation.R"
${R_PATH}/R --vanilla --slave --args ${MUTNUM} ${SAMPLENUM} ${ALPHA} ${GAMMA} ${OUTPUTFILE} < repeatSimulation.R
