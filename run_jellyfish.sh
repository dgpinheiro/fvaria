#!/bin/bash

for k in `seq 19 4 63`; do 
	jellyfish count -t 20 -C -m ${k} -s 10G -o ${k}mer_out --min-qual-char=? FVARIA_R1.fq FVARIA_R2.fq
	jellyfish histo -o ${k}mer_out.histo ${k}mer_out
done
