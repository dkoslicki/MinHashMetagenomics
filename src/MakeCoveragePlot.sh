#!/bin/bash
# This script will demonstrate how to call SNAP and samtools to generate the coverage figure
windowSize=10000
snap-aligner index PRJNA274798.fa . -s 20 -large
snap-aligner paired . 4539585.3.sorted.r1.fastq 4539585.3.sorted.r2.fastq -hp -mrl 40 -xf 1.2 -d 28 -o -sam out.sam
samtools sort aligned.sam > aligned.sorted.sam
samtools depth -a --reference PRJNA274798.fa out.sorted.sam | python get_coverage.py $windowSize /dev/fd/0 /dev/fd/1 > coverage_${windowSize}.txt
python CoveragePlot.py -i ../data/SNAP/coverage_${windowSize}.txt -o CoveragePlot.png