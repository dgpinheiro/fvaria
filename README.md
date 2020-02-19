# *Frieseomelitta varia* genome assembly

## Illumina HiSeq2500 sequencing:

The DNA library preparation and sequenging steps were made at [LaCTAD](https://www.lactad.unicamp.br/) sequencing facility.
Details about Illumina sequenging strategies and the resulting data are described:

* TruSeq Paired-end 
	* ~350 bp - average insert size
	* 2 * 101 bp
	* 1 Lane
	* 1 Sample
	
	| Name   | Lane | Index  | Sample   | #Reads          | % Bases >= Q30 |
	| ------ | ---- | ------ | -------- | --------------- | -------------- |
	| DNAPE1 | L006 | CAGATC | Drone    | 2 * 120,948,424 | 64.85          |
	| DNAPE2 | L001 | ACAGTG | Drone    | 2 * 130,859,645 | 94.53          |

* Nextera Mate-pair
	* ~3 Kbp - average fragment size
	* 2 * 101 bp
	* 1 Lane
	* 1 Sample
	
	| Name        | Lane | Index   | Sample   | #Reads          | % Bases >= Q30 |
	| ----------- | -----| ------- | -------- | --------------- | -------------- |
	| DNAMatePair | L007 | NoIndex | Drone    | 2 * 85,513,161  | 92.10          |

> The sequencing step was performed by LaCTAD and the above information was provided by MSc. Osvaldo Reis Júnior (LaCTAD).

## Genome size estimation
	
In order to execute genome size estimation based ok k-mer counting we need to execute Jellyfish ([run_jellyfish.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_jellyfish.sh)), but first we concatenated all the R1 files into one file (FVARIA_R1.fq) and all the R2 files into another file (FVARIA_R2.fq):

```bash=
cat ./raw/*_R1_*.fastq > FVARIA_R1.fq
cat ./raw/*_R2_*.fastq > FVARIA_R2.fq
```

```bash=
./run_jellyfish.sh
```

> Jellyfish was executed using a range of k-mer sizes, from 19 to 63 with an increment of 4. The Jellyfish histogram of k-mer counts was also created with [run_jellyfish.sh]().

* 1st method

	Genome size was estimated from all fastq reads by k-mer counting using [Jellyfish](https://github.com/gmarcais/Jellyfish) v. 2.2.0. 
The genome size estimation was based on [GenomeScope](https://github.com/schatzlab/genomescope) analysis. 

```bash=
for i in `ls *.histo`; do k=`basename ${i} .histo | sed 's/mer_out//'`; genomescope.R ${i} ${k} 101 ./GenomeScope${k}; done
```
```bash=
cat ./GenomeScope*/summary.txt | grep 'Genome Haploid Length' | perl -F"\s{2,}" -lane 'INIT { $min=0; $max=0; $n=0; } $n++; $F[1]=~s/\D+//g; $F[2]=~s/\D+//g; $min+=$F[1]; $max+=$F[2]; END { print $F[0],"\t",sprintf("%.2f",($min/$n)/1000000)," Mbp","\t",sprintf("%.2f", ($max/$n)/1000000)," Mbp"; }'
```

Output:
<pre>
Genome Haploid Length   345.77 Mbp      345.90 Mbp
</pre>

* 2st method

	Genome size was also estimated from all fastq reads by k-mer counting using [Jellyfish](https://github.com/gmarcais/Jellyfish) v. 2.2.0. 
The genome size estimation was based on [Genome Size Estimation Tutorial](https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/) of Computational Biology Core of Institute for Systems Biology at University of Connecticut. We create an R script ([genomeSize.R](https://github.com/dgpinheiro/fvaria/blob/master/genomeSize.R)) to automate the steps of this tutorial, including the identification and trimming the k-mers with frequencies less than the frequency of the first valley identified in the histogram of k-mer counts obtained with Jellyfish. The peak and the first valley was identified using [findPeaks](https://github.com/stas-g/findPeaks) function.
	The [genomeSize.R](https://github.com/dgpinheiro/fvaria/blob/master/genomeSize.R) was executed for a range o k-mer sizes and the average was calculated:

```bash=
perl -lane 'INIT {my $sum=0; my $n=0;} $n++; $sum+=$_; END { print "Genome size: ".sprintf("%.2f",($sum/$n))." Mbp"; }' <(for i in `ls *.histo`; do  ./genomeSize.R --in=./${i} | grep 'Genome size' | sed 's/Genome size: //' | sed 's/ Mb//'; done)
```

Output:
<pre>
Genome size: 389.51 Mbp
</pre>

* 3nd method

	Genome size was also estimated from all fastq reads by k-mer counting using [KMC](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) v. 3.0.0. To do this we set up kmc with the following parameters: to consider kmers of size 25 (-k25); to discard kmers with frequencies less than 5 (-ci5); to use 200 bins (-n200); and to process all using 20 threads (-t20). The file named "files.lst" contains the path to all fastq files. The number of unique counted k-mers was divided by 1,000,000 to obtain the estimated genome size in Mbp.

```bash=
kmc -k25 -ci5 -n200 -t20 @files.lst result /tmp
```

Output:
<pre>
Estimated Genome size (kmc): 381.41 Mbp
</pre>

## Pre-processing steps

The software [Fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to perform the data quality control.
Then the paired-end reads were analyzed with software [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
to identification of adapter sequences and quality filtering.

The mate-pair reads were analyzed using [NxTrim](https://github.com/sequencing/NxTrim) to
 to discard low quality reads and categorise reads according to the orientation implied by the adapter location
. Thus, it exploits the particular properties of Illumina Nextera Mate Pair data and evaluates the read on both sides          
of the adapter sites.
The NxTrim builds 'virtual libraries' of mate pairs, paired-end reads and single-ended reads. NxTrim
also trims off adapter read‐through.

* Post processing steps:

	| Protocol   | #Reads      |
	| ---------- | ----------- |
	| Paired-end | 234,357,438 |
	| Mate-pair  |   9,617,088 |

> The pre-processing step was performed by LaCTAD and the above information was provided by MSc. Osvaldo Reis Júnior (LaCTAD).

## Initial Genomic Assembly 

The assembly was performed using software [SPAdes](http://bioinf.spbau.ru/spades) v. 3.9.0
using error correction (BayesHammer module) and using 13 k-mer lengths (33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77 and 81).

<pre>
#Scaffolds: 9.755
Assembly size: 280.951.552 bp
Assembly size (without N's): 280.751.127
Average size of scaffolds: 28.800
N50: 81.551 (980 contigs)
</pre>

> The assembly step was performed by LaCTAD and the above information was provided by MSc. Osvaldo Reis Júnior (LaCTAD).

## Scaffolding Genomic Assembly

First, the raw mate-pairs reads were processed using [NxTrim](https://github.com/sequencing/NxTrim) according to the script [run_nxtrim.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_nxtrim.sh):

```bash=
./run_nxtrim.sh
```
> To execute the script [run_nxtrim.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_nxtrim.sh) we need the *./raw* directory containing the fastq files with *DNAPE* and *DNAMatePair* prefixes.

The software [BESST](https://github.com/ksahlin/BESST) v. 2.2.4 was used to increase scaffolding of the genomic assembly.
In order to do this, we first aligned all the processed reads to the assembled genome using [HISAT2](http://daehwankimlab.github.io/hisat2/) v. 2.0.5, and then we sorted and indexed the alignments using [samtools](http://samtools.sourceforge.net/) v. 1.3.1.

```bash=
cd ./ref

hisat2-build genome genome.fa

cd ../
```

> To execute the script [run_hisat2.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_hisat2.sh), first we need to create an index for the the fasta file with genome assembly (genome.fa) 
using software hisat2-build.

```bash=
./run_hisat2.sh
```

> The script [run_hisat2.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_hisat2.sh) must be executed in a directory containing the following input directories:
> - "*./in*" containing fastq files with *DNAPE* and *DNAMatePair* prefixes processed with [run_nxtrim.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_nxtrim.sh);
> - "*./ref*" containing hisat2 index files (\*.ht2);

So, we executed BESST pipeline with sorted and indexed alignments:

```bash=
mkdir ./out/

runBESST \
        -q \
        --orientation rf fr \
        --separate_repeats \
        --min_mapq 30 \
        -c ./ref/genome.fa \
	-g \
	-f \
        ./out/DNAMP_sorted.bam ./out/DNAPE_sorted.bam
```

> The final assembly was named *Fvar-1.1.fa*

## Gene prediction

The gene prediction was made with [Maker2](https://www.yandell-lab.org/software/maker.html) v. 3.00 using the *F. varia* assembled genome (*Fvar-1.1.fa*). In addition to the default gene prediction protocol in Maker2, we added the following data to be used:
* EST Evidence
	* *transcriptXgenome_filtered_sorted.gff* - *F. varia* assembled transcriptome (contigs) mapped on *F. varia* genome assembly (*Fvar-1.1.fa*) using [UCSC Genome Browser](https://genome-store.ucsc.edu/) Tools (blat, pslCDnaFilter, pslToBed, bedToGenePred, genePredToGtf) and [GBrowse](http://gmod.org/wiki/GBrowse) Tools (gtf2gff3);
	* *ests.fa* - 20,064 Expressed Sequence Tags (ESTs) of *F. varia* obtained from NCBI GenBank;
	* *altests.fa* - 169,511 Expressed Sequence Tags (ESTs) of *Apis mellifera* obtained from NCBI GenBank;
* Protein Homology Evidence
	* *GCF_000002195.4_Amel_4.5_protein.faa* - Proteins of *Apis mellifera* obtained from NCBI Genome;
Moreover, the following gene prediction softwares in Maker2 was also configured:
	* SNAP - configured with A.mellifera.hmm;
	* Augustus - configured with fly model;

> Full details can be obtained by consulting the Maker2 configuration files: maker_bopts.ctl  maker_evm.ctl  maker_exe.ctl  maker_opts.ctl
	
