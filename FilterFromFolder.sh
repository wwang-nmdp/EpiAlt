#!/bin/bash
FILES=/home/cdsw/MDS_bismark_cov_files/*.cov.gz

for f in $FILES
do
  echo "Removing low counts from $f files......."

python /home/cdsw/filterLowCounts.py $f "${f%.cov.gz}".filter10.txt 10

done
