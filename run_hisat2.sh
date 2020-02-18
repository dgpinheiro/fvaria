#!/bin/bash

if [ ! -d "./in/" ]; then
	echo "Error: Missing input directory: ./in containing .fastq files with DNAPE* and DNAMatePair* prefixes"
	exit
fi

dnape=()
for r1 in in/DNAPE*_R1.fastq.gz; do 
	bn=`basename ${r1} _R1.fastq.gz`
	r2=`echo ${r1} | sed 's/_R1./_R2./'`
	echo "Aligning ${bn} ..."
	dnape=(${dnape[@]} "out/${bn}.bam")
	time hisat2 --max-seeds 1000 --no-temp-splicesite --no-spliced-alignment --minins 100 --maxins 600 --end-to-end --threads 22 --fr -k 1 -x ./ref/genome -1 ${r1} -2 ${r2} 2> out/${bn}.log.txt | samtools view -b -S -o ./out/${bn}.bam -
done

samtools merge -f out/DNAPE.bam ${dnape[*]}

dnamp=()
for r1 in in/DNAMatePair_*_R1.mp.fastq.gz; do 
	bn=`basename ${r1} _R1.mp.fastq.gz`
	r2=`echo ${r1} | sed 's/_R1./_R2./'`
	echo "Aligning ${bn} ..."
	dnamp=(${dnamp[@]} "out/${bn}.bam")
	time hisat2 --max-seeds 1000 --no-temp-splicesite --no-spliced-alignment --minins 1000 --maxins 15000 --end-to-end --threads 22 --rf -k 1 -x ./ref/genome -1 ${r1} -2 ${r2} 2> out/${bn}.log.txt | samtools view -b -S -o ./out/${bn}.bam -
done

samtools merge -f ./out/DNAMP.bam ${dnamp[*]}

samtools sort ./out/DNAMP.bam -o ./out/DNAMP_sorted.bam
samtools index ./out/DNAMP_sorted.bam

