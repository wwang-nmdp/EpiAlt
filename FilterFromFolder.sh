#!/bin/bash
FILES=/home/wwang/MCW/bismark_out/*.cov.gz

for f in $FILES
do
  echo "Removing low counts from $f files......."

python /home/wwang/filterLowCounts.py $f "${f%.cov.gz}".filter10.txt 10

done
