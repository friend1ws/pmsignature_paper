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



echo -n > ../result/MPFormat/merge.mp.txt
while read type 
do 
    echo "zcat ../data/${type}_clean_somatic_mutations_for_signature_analysis.txt.gz | perl ../../common_script/convertToMPFormat.pl - >> ../result/MPFormat/merge.mp.txt"
    zcat ../data/${type}_clean_somatic_mutations_for_signature_analysis.txt.gz | perl ../../common_script/convertToMPFormat.pl - >> ../result/MPFormat/merge.mp.txt
done < ../data/cancer_types.txt

echo "gzip ../result/MPFormat/merge.mp.txt > ../result/MPFormat/merge.mp.txt.gz"
gzip ../result/MPFormat/merge.mp.txt > ../result/MPFormat/merge.mp.txt.gz

rm -rf ../result/MPFormat/merge.mp.txt

if [ ! -d log ]
then
    mkdir log
fi


SCRIPTDIR=`pwd`
cd ../result
RESULTDIR=`pwd`
cd ../data
DATADIR=`pwd`
cd ../../common_script


type=merge
for K in 5 10 15 20 25  
do

    echo "qsub -l s_vmem=8G,mem_req=8 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt.gz ${RESULTDIR}/Param_ind5/${type}.${K}.Rdata ${K} FALSE 1"
    qsub -l s_vmem=8G,mem_req=8 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt.gz ${RESULTDIR}/Param_ind5/${type}.${K}.Rdata ${K} FALSE 1

    echo "qsub -l s_vmem=8G,mem_req=8 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt.gz ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 1"
    qsub -l s_vmem=8G,mem_req=8 -e ${SCRIPTDIR}/log -o ${SCRIPTDIR}/log perform_pmsignature.sh ${RESULTDIR}/MPFormat/${type}.mp.txt.gz ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 1

done

