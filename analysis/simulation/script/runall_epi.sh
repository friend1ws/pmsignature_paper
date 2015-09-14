#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

if [ ! -d ../result ]
then
    mkdir ../result
fi

if [ ! -d ../result/cosineDist_epi ]
then
    mkdir ../result/cosineDist_epi
fi

if [ ! -d log ]
then
    mkdir log
fi

SCRIPTDIR=`pwd`
cd ../result
RESULTDIR=`pwd`
cd ../script/subscript
SUBSCRIPTDIR=`pwd`


for FEATNUM in 0 1 2 3 4 5; do
    # for SAMPLENUM in 10 25 50 100; do
    for SAMPLENUM in 10; do
        # for PARAM_ALPHA in 0.5 1 2; do
	for PARAM_ALPHA in 1; do
            # for PARAM_GAMMA in 0.5 1 2; do
	    for PARAM_GAMMA in 1; do

                OUTPUTFILE=${RESULTDIR}/cosineDist_epi/${FEATNUM}_${SAMPLENUM}_${PARAM_ALPHA}_${PARAM_GAMMA}.txt

                # echo "qsub -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log repeatSimulation_epi.sh ${FEATNUM} ${SAMPLENUM} ${PARAM_ALPHA} ${PARAM_GAMMA} ${OUTPUTFILE}"
                # qsub -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log repeatSimulation_epi.sh ${FEATNUM} ${SAMPLENUM} ${PARAM_ALPHA} ${PARAM_GAMMA} ${OUTPUTFILE}
	
		echo "sh repeatSimulation_epi.sh ${FEATNUM} ${SAMPLENUM} ${PARAM_ALPHA} ${PARAM_GAMMA} ${OUTPUTFILE}"
		sh repeatSimulation_epi.sh ${FEATNUM} ${SAMPLENUM} ${PARAM_ALPHA} ${PARAM_GAMMA} ${OUTPUTFILE}

            done;
        done;
    done;
done

