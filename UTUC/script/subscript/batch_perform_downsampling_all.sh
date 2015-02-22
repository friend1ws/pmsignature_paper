#! /bin/sh
#$ -S /bin/sh
#$ -cwd

source ../lib/config.sh

GROUP=$1
TYPE=$2

INPUT=${HOMEDIR}/data/UTUC/input/${GROUP}.${TYPE}.txt
OUTPUTDIR=${HOMEDIR}/data/UTUC/output/${GROUP}_${TYPE}

if [ ! -f ${INPUT} ]; then
    echo "${INPUT} does not exist."
    exit
fi

if [ ! -d ${OUTPUTDIR} ]; then
    mkdir -p ${OUTPUTDIR}
fi

for SIGNUM in `seq 2 7`; do
    if [ ! -f ${HOMEDIR}/data/output/${GROUP}_${TYPE}/${SIGNUM}_1_1.Rdata ]; then
        echo "qsub perform_downsampling.sh ${GROUP} ${TYPE} 1 ${SIGNUM} 1"
        qsub perform_downsampling.sh ${GROUP} ${TYPE} 1 ${SIGNUM} 1 
    fi;
done;


