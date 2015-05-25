#! /bin/sh
#$ -S /bin/sh
#$ -cwd


BEDTOOLS_PATH=/home/yshira/bin/bedtools2-2.22.0
GENOME_REF=/home/yshira/common/ref/hg19_all/hg19.all.fasta
PATH=${BEDTOOLS_PATH}/bin:$PATH

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/common_database/refGene.txt.gz
gunzip refGene.txt.gz

perl listRefGene.pl > ../../common_data/background/exon.bed
echo "bedtools random -l 5 -n 50000000 -g ${BEDTOOLS_PATH}/genomes/human.hg19.genome | sort -k1,1 -k2,2n > ../../common_common_data/background/rand.sorted.ed"
bedtools random -l 5 -n 50000000 -g ${BEDTOOLS_PATH}/genomes/human.hg19.genome | sort -k1,1 -k2,2n > ../../common_common_data/background/rand.sorted.bed

echo "sort -k1,1 -k2,2n ../../common_common_data/background/exon.bed > ../../common_common_data/background/exon.sorted.bed"
sort -k1,1 -k2,2n ../../common_common_data/background/exon.bed > ../../common_common_data/background/exon.sorted.bed

echo "bedtools intersect -a ../../common_common_data/background/rand.sorted.bed -b ../../common_common_data/background/exon.sorted.bed -s -wa -sorted | uniq > ../../common_common_data/background/rand.exon.bed"
bedtools intersect -a ../../common_common_data/background/rand.sorted.bed -b ../../common_common_data/background/exon.sorted.bed -s -wa -sorted | uniq > ../../common_common_data/background/rand.exon.bed

echo "bedtools getfasta -fi ${GENOME_REF} -bed ../../common_common_data/background/rand.exon.bed -fo ../../common_common_data/background/rand.exon.fasta"
bedtools getfasta -fi ${GENOME_REF} -bed ../../common_common_data/background/rand.exon.bed -fo ../../common_common_data/background/rand.exon.fasta 

echo "perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 0 > ../../common_data/background/bg.ind3.txt"
perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 0 > ../../common_data/background/bg.ind3.txt

echo "perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 1 > ../../common_data/background/bg.ind3_dir.txt"
perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 1 > ../../common_data/background/bg.ind3_dir.txt

echo "perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 0 > ../../common_data/background/bg.ind5.txt"
perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 0 > ../../common_data/background/bg.ind5.txt

echo "perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 1 > ../../common_data/background/bg.ind5_dir.txt"
perl proc_bg_ind.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 1 > ../../common_data/background/bg.ind5_dir.txt


echo "perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 0 > ../../common_data/background/bg.full3.txt"
perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 0 > ../../common_data/background/bg.full3.txt

echo "perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 1 > ../../common_data/background/bg.full3_dir.txt"
perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 3 1 > ../../common_data/background/bg.full3_dir.txt

echo "perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 0 > ../../common_data/background/bg.full5.txt"
perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 0 > ../../common_data/background/bg.full5.txt

echo "perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 1 > ../../common_data/background/bg.full5_dir.txt"
perl proc_bg_full.pl ../../common_data/background/rand.exon.fasta ../../common_data/background/rand.exon.bed 5 1 > ../../common_data/background/bg.full5_dir.txt


