#! /bin/sh
#$ -S /bin/sh
#$ -cwd

export R_LIBS=/home/yshira/.R
export R_PATH=/home/yshira/local/bin

if [ ! -d ../result/MFVFormat ]
then
    mkdir ../result/MFVFormat
fi

if [ ! -d ../result/Param_ind ]
then
    mkdir ../result/Param_ind
fi

if [ ! -d ../result/Param_full ]
then
    mkdir ../result/Param_full
fi


# prepare data
echo "perl subscript/proc.hoang.pl ../data/3006200s.txt > ../result/MFVFormat/Hoang_MFVF.temp.txt"
perl subscript/proc.hoang.pl ../data/3006200s.txt > ../result/MFVFormat/Hoang_MFVF.temp.txt

echo "perl ../../common_script/proc_ind.pl ../result/MFVFormat/Hoang_MFVF.temp.txt 5 1 > ../result/MFVFormat/Hoang_MFVF.ind.txt"
perl ../../common_script/proc_ind.pl ../result/MFVFormat/Hoang_MFVF.temp.txt 5 1 > ../result/MFVFormat/Hoang_MFVF.ind.txt 

echo "perl ../../common_script/proc_full.pl ../result/MFVFormat/Hoang_MFVF.temp.txt 5 1 > ../result/MFVFormat/Hoang_MFVF.full.txt"
perl ../../common_script/proc_full.pl ../result/MFVFormat/Hoang_MFVF.temp.txt 5 1 > ../result/MFVFormat/Hoang_MFVF.full.txt


# :<<_COMMENT_OUT_


if [ ! -d log ]
then
    mkdir log
fi

SCRIPTDIR=`pwd`
cd ../result
RESULTDIR=`pwd`
cd ../script/subscript
SUBSCRIPTDIR=`pwd`

SIGNUM=3
INPUTFILE=${RESULTDIR}/MFVFormat/Hoang_MFVF.ind.txt

for TYPE in ind full; do
    for RATIO in 0.01 0.025 0.05 0.1 0.25 0.5; do

        INPUTFILE=${RESULTDIR}/MFVFormat/Hoang_MFVF.${TYPE}.txt
        OUTPUTDIR=${RESULTDIR}/Param_${TYPE}
        echo "qsub -l s_vmem=2G,mem_req=2 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_downsampling.sh ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTDIR} 100 ${TYPE}"
        qsub -l s_vmem=2G,mem_req=2 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_downsampling.sh ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTDIR} 100 ${TYPE}

    done
done



TYPE=ind
RATIO=1
for SIGNUM in `seq 2 6`; do

    INPUTFILE=${RESULTDIR}/MFVFormat/Hoang_MFVF.${TYPE}.txt
    OUTPUTDIR=${RESULTDIR}/Param_${TYPE}
    qsub -l s_vmem=2G,mem_req=2 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_downsampling.sh ${INPUTFILE} ${RATIO} ${SIGNUM} ${OUTPUTDIR} 1 ${TYPE}

done

# _COMMENT_OUT_


