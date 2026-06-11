#!/bin/bash

for f in $(ls *fastq.gz)
do
	id=$(echo $f | cut -f 1 -d ".")
	NanoPlot -t 12 --tsv_stats --N50 --fastq $f -o ${id}_raw
	zcat $f	| NanoFilt --q 8 -l 1200 --maxlength 1800 | gzip > ${id}_filt.fastq.gz
	NanoPlot -t 12 --tsv_stats --N50 --fastq ${id}_filt.fastq.gz -o ${id}_filt
done

