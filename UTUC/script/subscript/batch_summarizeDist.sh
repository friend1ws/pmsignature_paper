#! /bin/sh
#$ -S /bin/sh
#$ -cwd

GROUP=$1
TYPE=$2
SIGNUM=$3

source ../lib/config.sh


for RATIO in 0.01 0.025 0.05 0.1 0.25 0.5; do

        INPUTDIR=${HOMEDIR}/data/UTUC/output/${GROUP}_${TYPE}/${SIGNUM}_${RATIO}
        TRUEFILE=${HOMEDIR}/data/UTUC/output/${GROUP}_${TYPE}/${SIGNUM}_1_1.Rdata
        OUTPUTDIR=${HOMEDIR}/data/UTUC/output/summary/${GROUP}_${TYPE}_${SIGNUM}_${RATIO}.txt

        echo "${R_PATH}/R --vanilla --slave --args ${INPUTDIR} ${TRUEFILE} ${OUTPUTDIR} ${TYPE} < summarizeDist.R"
        ${R_PATH}/R --vanilla --slave --args ${INPUTDIR} ${TRUEFILE} ${OUTPUTDIR} ${TYPE} < summarizeDist.R
done

