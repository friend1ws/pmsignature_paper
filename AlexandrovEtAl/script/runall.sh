#! /bin/sh

if [ ! -d ../result/MPFormat ]
then
    mkdir ../result/MPFormat
fi

if [ ! -d ../result/Param_ind5 ]
then
    mkdir ../result/Param_ind5
fi

if [ ! -d ../result/Param_ind5_dir ]
then
    mkdir ../result/Param_ind5_dir
fi


:<<_COMMENT_OUT_

while read type 
do 
    echo "perl ../../common_script/convertToMPFormat.pl ../data/AlexandrovEtAl_raw/${type}_clean_somatic_mutations_for_signature_analysis.txt > ../result/MPFormat/${type}.mp.txt"
    perl ../../common_script/convertToMPFormat.pl ../data/AlexandrovEtAl_raw/${type}_clean_somatic_mutations_for_signature_analysis.txt > ../result/MPFormat/${type}.mp.txt
done < ../data/AlexandrovEtAl_raw/cancer_types.txt

_COMMENT_OUT_


cd ../result
RESULTDIR=`pwd`
cd ../data
DATADIR=`pwd`
cd ../../common_script
SCRIPTDIR=`pwd`

while read type
do

    for K in `seq 2 6`
    do
        echo "qsub -l s_vmem=2G,mem_req=2 perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5/${type}.${K}.Rdata ${K} FALSE 10"
        qsub -l s_vmem=2G,mem_req=2 perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5/${type}.${K}.Rdata ${K} FALSE 10

        if [ ${type} == Lung-Adeno ]
        then 
            echo "qsub -l s_vmem=4G,mem_req=4 perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 10"
            # qsub -l s_vmem=4G,mem_req=4 perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 10
        else
            echo "qsub -l s_vmem=2G,mem_req=2 perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 10"
            # qsub -l s_vmem=2G,mem_req=2 perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 10
        fi

    done

done < ${DATADIR}/AlexandrovEtAl_raw/cancer_types.txt

