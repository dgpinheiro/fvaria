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

## Genome estimates 

The genome heterozygosity, repeat content, and size was evaluated from sequencing reads using a kmer-based statistical approach using [GenomeScope](https://github.com/schatzlab/genomescope) software v. 1.0.

In order to perform the estimates we need to execute [Jellyfish](https://github.com/gmarcais/Jellyfish) v. 2.2.0 ([run_jellyfish.sh](https://github.com/dgpinheiro/fvaria/blob/master/run_jellyfish.sh)) for counting of k-mers in DNA. So, we first concatenated all the R1 files into one file (FVARIA_R1.fq) and all the R2 files into another file (FVARIA_R2.fq):

```bash=
cat ./raw/*_R1_*.fastq > FVARIA_R1.fq
cat ./raw/*_R2_*.fastq > FVARIA_R2.fq
```
	
```bash=
./run_jellyfish.sh
```

> Jellyfish was executed using a range of k-mer sizes, from 19 to 63 with an increment of 4. The Jellyfish histogram of k-mer counts was also created with [run_jellyfish.sh]().

Execution of GenomeScope for all Jellyfish histogram files:

```bash=
for i in `ls *.histo`; do k=`basename ${i} .histo | sed 's/mer_out//'`; genomescope.R ${i} ${k} 101 ./GenomeScope${k}; done
```

* Genome size

```bash=
cat ./GenomeScope*/summary.txt | grep 'Genome Haploid Length' | perl -F"\s{2,}" -lane 'INIT { $min=0; $max=0; $n=0; } $n++; $F[1]=~s/\D+//g; $F[2]=~s/\D+//g; $min+=$F[1]; $max+=$F[2]; END { print $F[0],"\t",sprintf("%.2f",($min/$n)/1000000)," Mbp","\t",sprintf("%.2f", ($max/$n)/1000000)," Mbp"; }'
```
Output (Min./Max.):
<pre>
Genome Haploid Length   345.77 Mbp      345.90 Mbp
</pre>

* Heterozygosity

The average of minimum and maximum heterozygosity estimates using distinct k-mer sizes:

```bash=
cat ./GenomeScope*/summary.txt | grep 'Heterozygosity' | perl -F"\s{2,}" -lane 'INIT { $min=0; $max=0; $n=0; } $n++; $F[1]=~s/\%//g; $F[2]=~s/\%//g; $min+=$F[1]; $max+=$F[2]; END { print $F[0],"\t",sprintf("%.4f",($min/$n))," %","\t",sprintf("%.4f", ($max/$n))," %"; }'
```

Output (Min./Max.):
<pre>
Heterozygosity  0.1258 %  0.1277 %
</pre>

* Genome Repeat Length

The average of minimum and maximum genome repeat length estimates using distinct k-mer sizes:

```bash=
cat ./GenomeScope*/summary.txt | grep 'Genome Repeat Length' | perl -F"\s{2,}" -lane 'INIT { $min=0; $max=0; $n=0; } $n++; $F[1]=~s/\D+//g; $F[2]=~s/\D+//g; $min+=$F[1]; $max+=$F[2]; END { print $F[0],"\t",sprintf("%.2f",($min/$n)/1000000)," Mbp","\t",sprintf("%.2f", ($max/$n)/1000000)," Mbp"; }'
```

Output (Min./Max.):
<pre>
Genome Repeat Length    76.05 Mbp       76.08 Mbp
</pre>

* Model Fit

The average of minimum and maximum model fit estimates using distinct k-mer sizes:

```bash=
cat ./GenomeScope*/summary.txt | grep 'Model Fit' | perl -F"\s{2,}" -lane 'INIT { $min=0; $max=0; $n=0; } $n++; $F[1]=~s/\%//g; $F[2]=~s/\%//g; $min+=$F[1]; $max+=$F[2]; END { print $F[0],"\t",sprintf("%.4f",($min/$n))," %","\t",sprintf("%.4f", ($max/$n))," %"; }'
```

Output (Min./Max.):
<pre>
Model Fit       97.1158 %       99.3401 %
</pre>

* Read Error Rate

The average of minimum and maximum read error rate estimates using distinct k-mer sizes:

```bash=
cat ./GenomeScope*/summary.txt | grep 'Read Error Rate' | perl -F"\s{2,}" -lane 'INIT { $min=0; $max=0; $n=0; } $n++; $F[1]=~s/\%//g; $F[2]=~s/\%//g; $min+=$F[1]; $max+=$F[2]; END { print $F[0],"\t",sprintf("%.8f",($min/$n))," %","\t",sprintf("%.8f", ($max/$n))," %"; }'
```

Output (Min./Max.):
<pre>
Read Error Rate 0.12025125 %    0.12025125 %
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
#Scaffolds: 9,755
Assembly size: 280,951,552 bp
Average size of scaffolds: 28,800 bp
N50: 81,551 (980 contigs)
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

> The final assembly was named *Fvar-1.1.fa* and then copied to *Fvar-1.2.fa*, i.e. these fasta files are the same. The genome fasta file was compressed in [Fvar-1.2.fa.gz](https://github.com/dgpinheiro/fvaria/blob/master/data/Fvar-1.2.fa.gz).

<pre>
#Scaffolds: 2,173
Assembly size: 275,421,029 bp
Average size of scaffolds: 126,747 bp
N50: 467,007 (176 contigs)
</pre>


## Gene prediction

The gene prediction was made with [Maker2](https://www.yandell-lab.org/software/maker.html) v. 3.00 using the *F. varia* assembled genome (*Fvar-1.1.fa*).
First, we used Maker2 with the default gene prediction protocol with *Ab-initio* gene prediction softwares and some EST Evidences.

* EST Evidence
	* *ests.fa* - 20,064 Expressed Sequence Tags (ESTs) of *F. varia* obtained from NCBI GenBank;
	* *altests.fa* - 169,511 Expressed Sequence Tags (ESTs) of *Apis mellifera* obtained from NCBI GenBank;
* Protein Homology Evidence
	* *GCF_000002195.4_Amel_4.5_protein.faa* - Proteins of *Apis mellifera* obtained from NCBI Genome;
Moreover, the following *Ab-initio* gene prediction softwares in Maker2 was also configured:
	* [SNAP](https://github.com/KorfLab/SNAP) v. 2006-07-28 - configured with A.mellifera.hmm;
	* [Augustus](http://bioinf.uni-greifswald.de/augustus/) v. 3.3 - configured with fly model;
	* [genemark](http://exon.gatech.edu/GeneMark/) v. 3.52 - eukaryotic executable;

> Full details can be obtained by consulting the Maker2 configuration files ([maker2_1](https://github.com/dgpinheiro/fvaria/tree/master/data/maker2_1)): maker_bopts.ctl  maker_evm.ctl  maker_exe.ctl  maker_opts.ctl
 
So, we executed Maker2 again with another gene evidences. In addition to the previous configuration, we added another EST evidence based on transcriptome assembly together with the previous assembly:

* Re-annotation Using MAKER Derived GFF3
	* *Fvar-1.1.maker.output/Fvar-1.1.all.gff* - previous MAKER Derived GFF3;
* EST evidence
	* *transcriptXgenome_filtered_sorted.gff* - *F. varia* assembled transcriptome (contigs) mapped on *F. varia* genome assembly (*Fvar-1.1.fa*) using [UCSC Genome Browser](https://genome-store.ucsc.edu/) Tools (blat, pslCDnaFilter, pslToBed, bedToGenePred, genePredToGtf) and [GBrowse](http://gmod.org/wiki/GBrowse) Tools (gtf2gff3);

> Full details can be obtained by consulting the Maker2 configuration files ([maker2_2](https://github.com/dgpinheiro/fvaria/tree/master/data/maker2_2)): maker_bopts.ctl  maker_evm.ctl  maker_exe.ctl  maker_opts.ctl

## Creating GFF/GTF

We create the gene feature file from the output of Maker2 (an output folder with a datastore index) using *gff3_merge* script to combine all the GFFs.

```bash=
gff3_merge -d Fvar-1.2.maker.output/Fvar-1.2_master_datastore_index.log
```

We filtered out all the complementary features contained in the generated GFF file.

```bash=
grep -v -P '\t(contig|repeatmasker|snap_masked|augustus_masked_match|augustus_masked|blastx|protein2genome|evm|blastn|est2genome|cdna2genome|repeatrunner|tblastx)\t' Fvar-1.2.all.gff | grep -v '^[#>]' | grep -v -P '^[ACGTN]+$' > Fvar-1.2.cleaned.gff
```

We also converted GFF to GTF file:

```bash=
cat Fvar-1.2.cleaned.gff | gff2gtf.pl > Fvar-1.2.cleaned.gtf
```

So, we renamed the features using the prefix **Fvar** using first the [rename_gtf.pl](https://github.com/dgpinheiro/fvaria/blob/master/rename_gtf.pl) script and then the [rename_gff.pl](https://github.com/dgpinheiro/fvaria/blob/master/rename_gff.pl) script. After renaming with *rename_gff.pl* the [checkGFFcoords.pl](https://github.com/dgpinheiro/fvaria/blob/master/checkGFFcoords.pl) script was executed to check *three_prime_UTR* and *five_prime_UTR* features and adjust the first/last exon coordinates.

```bash=
rename_gtf.pl -i Fvar-1.2.cleaned.gtf -o Fvar-1.2.gtf -p Fvar

./rename_gff.pl -i Fvar-1.2.cleaned.gff -o Fvar-1.2.tmp.gff
./checkGFFcoords.pl -i1 Fvar-1.2.all.gff -i2 Fvar-1.2.tmp > Fvar-1.2.gff

rm -f ./Fvar-1.2.tmp.gff
```

> All these steps were made for the first execution of Maker2, but only the last one was represented in the commands above.


## Annotation with eggNOG-mapper

A coding gene function annotation was performed with [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) v. 2.0.0 before submission to NCBI Genome Database. The protein sequences (*Fvar-1.2-proteins.fa*) were extracted from GFF (*Fvar-1.2.gff*) using [gffread](https://github.com/gpertea/gffread) software v. 0.11.7.

```bash=
mkdir -p em/

sed 's/\.$//' Fvar-1.2-proteins.fa > Fvar-1.2-proteins-cleaned.fa

perl -e 'use Bio::SeqIO; my $seqin = Bio::SeqIO->new(-file=>"Fvar-1.2-proteins-cleaned.fa", -format=>"FASTA"); my $seqout = Bio::SeqIO->new(-fh=>\*STDOUT,-format=>"FASTA",-width=>1000000); while (my $seq=$seqin->next_seq() ) { $seqout->write_seq($seq); } ' | split -l 20000 -a 3 -d - em/input_file.chunk_

for f in em/*.chunk_*; do
       emapper.py -m diamond --no_annot --no_file_comments --cpu 40 -i ${f} -o ${f};
done
cat em/*.chunk_*.emapper.seed_orthologs > em/input_file.emapper.seed_orthologs

emapper.py --annotate_hits_table em/input_file.emapper.seed_orthologs --no_file_comments -o em/output_file --cpu 40

```

## UCSC Genome Browser

We load *F. varia* genome and gene prediction coordinates to our local UCSC Genome Browser in [*F. varia*](http://kerr.fmrp.usp.br:88/cgi-bin/hgGateway?hgsid=5040&clade=insect&org=F.+varia&db=0) page.

## Assembly of Mitochondrial Genome



## NCBI Genome Database Submission

* Nuclear genome

We first adjusted the GFF to satisfy NCBI genome submission requirements, such as locus_tag attribute. We also added eggNOG-mapper annotation for the coding genes.

```bash=
./adj_gff.pl -i Fvar-1.2.gff -e ./output_file.emapper.annotations > Fvar-1.2-ncbi.gff
```

The GFF obtained with Maker2 was annotated with [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) v. 
The GFF (*Fvar-1.2-ncbi.gff*) was converted to Sequin data file (*Fvar-1.2.sqn*) and was sent to NCBI Genome database.

```bash=
table2asn -M n -J -c w -euk -t ./template.sbt -gaps-min 10 -l paired-ends -i ./Fvar-1.2.fa -f ./Fvar-1.2-ncbi.gff -o ./Fvar-1.2.sqn  -n "Frieseomelitta varia" -taxid 561572 -V b -T -j "[organism=Frieseomelitta varia] [topology=linear] [location=genomic] [moltype=DNA] [gcode=1] [sex=male] [country=Brazil] [dev-stage=pharate-adult] [Tech=wgs] [completedness=partial] [common=marmelada] [strand=double]"
```

* Mitochondrial genome

The TBL file (*Fvar-MT-1.2.1.tbl*) obtained with [Sequin](https://www.ncbi.nlm.nih.gov/Sequin/) was manually adjusted (for example, the names of gene products were adjusted to follow the pattern established by NCBI), then converted to Sequin file (*Fvar-MT-1.2.1.sqn*) and was also sent to NCBI Genome database.

```bash=
table2asn -C lbda -Z -M n -J -c 'wsdD' -euk -t ./template.sbt -gap-type scaffold -l paired-ends -i ./Fvar-MT-1.2.fa -f ./Fvar-MT-1.2.1.tbl -o ./Fvar-MT-1.2.1.sqn  -n "Frieseomelitta varia" -taxid 561572 -V b -T -j "[organism=Frieseomelitta varia] [topology=circular] [location=mitochondrion] [moltype=DNA] [mgcode=5] [gcode=1] [sex=male] [country=Brazil] [dev-stage=pharate-adult] [tech=wgs] [completedness=complete] [common=marmelada] [strand=double]"
```

