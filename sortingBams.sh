#!/bin/bash
FILES=/home/wwang/MCW/BAM/*.bam

for f in $FILES
do
  echo "Sorting $f ......."

samtools view -h $f > /home/wwang/MCW/SAMs/${f##*/}.sam && samtools sort /home/wwang/MCW/SAMs/${f##*/}.sam > /home/wwang/MCW/SAMs/${f##*/}.sorted.bam && samtools index /home/wwang/MCW/SAMs/${f##*/}.sorted.bam && rm /home/wwang/MCW/SAMs/${f##*/}.sam

done
