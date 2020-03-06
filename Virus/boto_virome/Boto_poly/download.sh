#!/bin/bash

#Script to download Botocudo files.
#In the file "samples2download.txt" the link to download for each sample is written
i=2
while read line; do
	wget -o log_${i}.txt $line
	i=$((i+1))
done <samples2download.txt

wait

for j in *.fastq.gz; do
	md5sum $j >> checksums.txt 
done
