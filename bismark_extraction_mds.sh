#!/bin/bash
FILES=/home/wwang/MCW/BAM/*.bam

for f in $FILES
do
  echo "extracting CpG $f ......."
perl /home/wwang/Tools/Bismark-0.20.0/bismark_methylation_extractor -p --no_overlap -o /home/wwang/MCW/bismark_out/ --merge_non_CpG --comprehensive --bedGraph --ignore 5 --counts --report $f && rm /home/wwang/MCW/bismark_out/Non_CpG_context_MCW2018*


done &&

FILES=/home/wwang/MCW/SAMs/*.bam

for f in $FILES
do
  echo "Calculate the mean read depth of $f ......."
samtools depth -a $f | awk '{c++;s+=$3}END{print s/c}'

done &&

for f in $FILES
do
  echo "Calculate the breadth coverage of $f ......."
samtools depth -a $f | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'

done &&

for f in $FILES
do
  echo "Proportion of the Reads that Mapped to the Reference of $f ......."
samtools flagstat $f | awk -F "[(|%]" 'NR == 3 {print $2}'

done
