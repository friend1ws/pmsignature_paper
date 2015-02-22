#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

INPUTFILE=$1
OUTPUTFILE=$2
SIGNUM=$3
TRDIRFLAG=$4
TRIALNUM=$5

echo "${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${OUTPUTFILE} ${SIGNUM} ${TRDIRFLAG} ${TRIALNUM} < perform_AlexandrovEtAl.R"
${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${OUTPUTFILE} ${SIGNUM} ${TRDIRFLAG} ${TRIALNUM} < perform_AlexandrovEtAl.R

