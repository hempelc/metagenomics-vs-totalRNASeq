#!/bin/bash

num_fil=$(find . -maxdepth 1 -type f | wc -l)
echo "Number of files: $num_fil"

failed=$(du *final.txt | grep -P "^4\t" | wc -l)
echo "Number of failed pipelines: $failed"

if [[ $failed != 0 ]]; then
  echo "Generate file containing failed pipelines: failed.txt"
  du *final.txt | grep -P "^4\t" | cut -f 2 | sed 's/\/.*$//g' > failed.txt
fi
