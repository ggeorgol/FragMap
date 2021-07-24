#!/bin/bash

####	FragMap version 2.0 Workflow to Map Sequencing Reads to Genomic Locations	 ####
####	 		©Grigorios Georgolopoulos 2017-2019		####
module load bcl2fastq2
module load jdk
module load bwa
module load samtools
module load bedops

FCPATH=$1
REF=$2
ADAPTERS=$3

cd "$FCPATH"

mkdir "$FCPATH"/Trimmed
TRIMPATH=$FCPATH/Trimmed

mkdir "$FCPATH"/fastq_files
echo $(date -R) "Generating fastq files"

bcl2fastq -R $FCPATH -o $FCPATH/fastq_files --adapter-stringency 0.75

cd $FCPATH/fastq_files

echo $(date -R) "Trimming adapter sequences"

for f in $(ls | grep R1 | grep -v Und)
do
	/home/ggeorgol/.local/bin/seqtk/seqtk trimfq -b 22 $f | gzip -c > $TRIMPATH/$f.fq.gz
done

for f in $(ls | grep R2 | grep -v Und)
do
	/home/ggeorgol/.local/bin/seqtk/seqtk trimfq -b 20 $f | gzip -c > $TRIMPATH/$f.fq.gz
done

cd "$TRIMPATH"

shopt -s nullglob

r1=(*R1*) # Read 1 files
r2=(*R2*) # Read 2 files
str=($(ls | grep R1 | cut -c1-4)) # Get the first 3 characters for output filenames

tLen=${#r1[@]}

for (( i=0; i<${tLen}; i++ ));
do
        echo $(date -R) "Processing ${str[$i]}"
        java -jar /home/ggeorgol/.local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 16 $TRIMPATH/${r1[$i]} $TRIMPATH/${r2[$i]} -baseout $TRIMPATH/${str[$i]} ILLUMINACLIP:${ADAPTERS}.fa:2:30:10:2:true
done

mkdir $TRIMPATH/out
OUTPATH=$TRIMPATH/out

for (( i=0; i<${tLen}; i++ ));
do
        echo $(date -R) "Aligning ${str[$i]}"
        bwa mem /net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts $TRIMPATH/${str[$i]}_1P $TRIMPATH/${str[$i]}_2P > $OUTPATH/${str[$i]}.sam
		samtools view -b $OUTPATH/${str[$i]}.sam | samtools sort - > $OUTPATH/${str[$i]}.bam
		samtools index $OUTPATH/${str[$i]}.bam
		echo $(date -R) "Converting ${str[$i]}.bam to ${str[$i]}.bed"
		samtools view -h -q 1 -f 64 $OUTPATH/${str[$i]}.bam | bam2bed - > $OUTPATH/${str[$i]}.bed
done

echo $(date -R) "Alignment completed"

for f in `ls $OUTPATH/*bed | grep GG` ; do ID=`basename $f | cut -f1 -d '.' | sed 's/_//'` ; IDnumber=`echo $ID | sed 's/[A-Z]\+//'` ; IDletter=`echo $ID | sed 's/[0-9]//g'` ; bedmap --delim "\t" --fraction-map 0.9 --echo --count $REF $f | awk -v IDletter=$IDletter -v IDnumber=$IDnumber 'BEGIN{OFS="\t"}{print $1,$2,$3,IDletter,IDnumber,$4}' >> $OUTPATH/output.txt ; done ; sort -k1,1 -k2n,2n -k3n,3n -k4,4 -k5n,5n $OUTPATH/output.txt | awk 'BEGIN{prevChr="xxx";prevBeg=-1;prevEnd=-1}{chr=$1;beg=$2;end=$3;if(chr != prevChr || beg != prevBeg || end != prevEnd){if(NR>1){printf("\n")}printf("%s\t%s\t%s\t%s",chr,beg,end,$6)}else{printf("\t%s",$6)}prevChr=chr;prevBeg=beg;prevEnd=end}END{printf("\n")}' > $OUTPATH/output.bed ; rm -f $OUTPATH/output.txt

echo $(date -R) "FragMap Completed Successfully"
