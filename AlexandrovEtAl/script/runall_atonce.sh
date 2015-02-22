#! /bin/sh

if [ ! -d MPFormat ]
then
    mkdir MPFormat
fi

if [ ! -d Param_ind5 ]
then
    mkdir Param_ind5
fi

if [ ! -d Param_ind5_dir ]
then
    mkdir Param_ind5_dir
fi



echo -n > MPFormat/all.mp.txt
while read type 
do 
    echo "perl ../script/command_AlexandrovEtAl/convertToMPFormat.pl ../data/AlexandrovEtAl_raw/${type}_clean_somatic_mutations_for_signature_analysis.txt >> MPFormat/all.mp.txt"
    perl ../script/command_AlexandrovEtAl/convertToMPFormat.pl ../data/AlexandrovEtAl_raw/${type}_clean_somatic_mutations_for_signature_analysis.txt >> MPFormat/all.mp.txt
done < ../data/AlexandrovEtAl_raw/cancer_types.txt


RESULTDIR=`pwd`
cd ../script/command_AlexandrovEtAl
SCRIPTDIR=`pwd`


type=all
for K in 5 10 15 20 25  
do

    echo "qsub -l s_vmem=8G,mem_req=8 perform_AlexandrovEtAl.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5/${type}.${K}.Rdata ${K} FALSE 1"
    qsub -l s_vmem=8G,mem_req=8 perform_AlexandrovEtAl.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5/${type}.${K}.Rdata ${K} FALSE 1

    echo "qsub -l s_vmem=8G,mem_req=8 perform_AlexandrovEtAl.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 1"
    qsub -l s_vmem=8G,mem_req=8 perform_AlexandrovEtAl.sh ${RESULTDIR}/MPFormat/${type}.mp.txt ${RESULTDIR}/Param_ind5_dir/${type}.${K}.Rdata ${K} TRUE 1

done

