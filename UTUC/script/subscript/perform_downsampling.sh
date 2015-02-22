#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

INPUTFILE=$1
RATIO=$2
SIGNUM=$3
OUTPUTDIR=$4
REPNUM=$5
TYPE=$6

echo "${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTDIR} ${REPNUM} ${TYPE} < downSampling.R"
${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTDIR} ${REPNUM} ${TYPE} < downSampling.R

