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

