# *Frieseomelitta varia* genome assembly


## Illumina HiSeq2500 sequencing:

The DNA library preparation and sequenging steps were made at [LaCTAD](https://www.lactad.unicamp.br/) sequencing facility.
Details about Illumina sequenging strategies and the resulting data are described:

* Paired-end
	* 2 * 101 bp
	* 1 Lane
	* 1 Sample

	| Lane | Sample   | #Reads      | % Bases >= Q30 |
	| ---- | -------- | ----------- | -------------- |
	| 1    | Amostra2 | 261.719.290 | 94,53          |

* Mate-pair
	* 2 * 101 bp
        * 1 Lane
        * 1 Sample

	| Lane | Sample   | #Reads      | % Bases >= Q30 |
        | -----| -------- | ----------- | -------------- |
	| 2    | Amostra1 | 171.026.322 | 90,00          |

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
	| Paired-end | 234.357.438 |
        | Mate-pair  |   9.617.088 |

## Assembly

The assembly was performed using software [SPAdes](http://bioinf.spbau.ru/spades) v. 3.9.0
using error correction (BayesHammer module) and using 13 k-mer lengths (33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77 and 81).

<pre>
#Scaffolds: 9.755
Assembly size: 280.951.552 bp
Assembly size (without N's): 280.751.127
Average size of scaffolds: 28.800
N50: 81.551 (980 contigs)
</pre>


