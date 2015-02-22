#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -e ../../log/ -o ../../log/

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

INPUTFILE=$1
OUTPUTFILE=$2
SIGNUM=$3
TRDIRFLAG=$4
TRIALNUM=$5

if [ ! -f ${INPUT} ]; then
    echo "${INPUT} does not exist."
    exit
fi

if [ ! -d ${OUTPUTDIR} ]; then
    mkdir -p ${OUTPUTDIR}
fi


  6 inputFile <- commandArgs()[5];
  7 ratio <- as.numeric(commandArgs()[6]);
  8 sigNum <- as.integer(commandArgs()[7]);
  9 outputFile <- commandArgs()[8];

for RATIO in 0.01 0.025 0.05 0.1 0.25 0.5; do
    for REPNUM in `seq 1 100`; do

        OUTPUTFILE=${OUTPUTDIR}/${SUGNUM}.${RATIO}.${REPNUM}.Rdata
        echo "${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTFILE} < downSampling.R"
        ${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTFILE} < downSampling.R

    done
done

GROUP=$1
TYPE=$2
SIGNUM=$3

source ../lib/config.sh

INPUT=${HOMEDIR}/data/UTUC/input/${GROUP}.${TYPE}.txt
OUTPUTDIR=${HOMEDIR}/data/UTUC/output/${GROUP}_${TYPE}

if [ ! -f ${INPUT} ]; then
    echo "${INPUT} does not exist."
    exit
fi

if [ ! -d ${OUTPUTDIR} ]; then
    mkdir -p ${OUTPUTDIR}
fi


for RATIO in 0.01 0.025 0.05 0.1 0.25 0.5; do
    for REPNUM in `seq 1 100`; do

        echo "qsub perform_downsampling.sh ${GROUP} ${TYPE} ${RATIO} ${SIGNUM} ${REPNUM}"
        qsub perform_downsampling.sh ${GROUP} ${TYPE} ${RATIO} ${SIGNUM} ${REPNUM} 

    done;
done

