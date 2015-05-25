#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

INPUTFILE=$1
OUTPUTFILE=$2
SIGNUM=$3
TYPE=$4

echo "${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${OUTPUTFILE} ${SIGNUM} ${TYPE} < perform_pmsignature_MFV.R"
${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${OUTPUTFILE} ${SIGNUM} ${TYPE} < perform_pmsignature_MFV.R

