#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

if [ ! -d ../result ]
then
    mkdir ../result
fi

if [ ! -d ../result/cosineDist ]
then
    mkdir ../result/cosineDist
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


for MUTNUM in 10 25 50 100 250 500 1000; do
    for SAMPLENUM in 10 25 50 100; do
        for PARAM_ALPHA in 0.5 1 2; do
            for PARAM_GAMMA in 0.5 1 2; do

                OUTPUTFILE=${RESULTDIR}/cosineDist/${MUTNUM}_${SAMPLENUM}_${PARAM_ALPHA}_${PARAM_GAMMA}.txt

                echo "qsub -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log repeatSimulation.sh ${MUTNUM} ${SAMPLENUM} ${PARAM_ALPHA} ${PARAM_GAMMA} ${OUTPUTFILE}"
                qsub -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log repeatSimulation.sh ${MUTNUM} ${SAMPLENUM} ${PARAM_ALPHA} ${PARAM_GAMMA} ${OUTPUTFILE}

            done;
        done;
    done;
done

