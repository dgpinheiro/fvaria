#!/bin/bash

if [ ! -d "./raw/" ]; then
	echo "Error: Missing input directory: ./in containing .fastq files with DNAPE* and DNAMatePair* prefixes"
	exit
fi

if [ ! -d "./in/" ]; then
	echo "Error: Missing input directory: ./in containing .fastq files with DNAPE* and DNAMatePair* prefixes"
	exit
fi

for r1 in ./raw/DNAMatePair*_R1*.fastq; do
	bn=`basename ${r1} .fastq | sed 's/_R1.*/_R1/'`;
	echo ${bn}; 
	cat ${r1} >> ./in/${bn}.fastq; 
done

echo "Concatenating DNAMatePair*.fastq files ..."

for r2 in ./raw/DNAMatePair*_R2*.fastq; do 
	bn=`basename ${r2} .fastq | sed 's/_R2.*/_R2/'`; 
	echo ${bn}; cat ${r2} >> ${bn}.fastq;
done

echo "Analysis with NxTrim ..."

nxtrim 	-1 ./in/DNAMatePair_NoIndex_L007_R1.fastq \
	-2 ./in/DNAMatePair_NoIndex_L007_R2.fastq \
	--rf \
	--separate \
	--aggressive \
	--similarity 0.85 \
	--minoverlap 8 \
	--minlength 21 \
	--output-prefix DNAMatePair_NoIndex_L007 \
	&> DNAMatePair_NoIndex_L007.log.txt

echo "Compressing Mate-pair files ..."

gzip ./in/DNAMatePair_NoIndex_L007_R1.fastq
gzip ./in/DNAMatePair_NoIndex_L007_R2.fastq

