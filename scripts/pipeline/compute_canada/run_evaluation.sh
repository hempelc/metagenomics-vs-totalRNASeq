#!/bin/bash

cluster=$1

num_dir=$(find . -maxdepth 1 -type d -iname "[0-9]*" | wc -l)
echo "Number of directories: $num_dir"

if [[ $num_dir != 512 ]]; then
  echo "Generate file containing missing pipelines: missing.txt"
  cp /home/hempelc/chrisnatjulia/scripts/pipeline/compute_canada/combinations_metagenomics_metatranscriptomics_pipeline_${cluster}.txt .
  mv combinations_metagenomics_metatranscriptomics_pipeline_${cluster}.txt all.txt
  touch done.txt
  for i in */; do
    echo $i >> done.txt
    sed -i 's/\///g' done.txt
  done
  cat done.txt all.txt | sort | uniq -u > missing.txt
  rm all.txt done.txt
fi

failed_started=$(du */ | grep -P "^4\t" | wc -l)
failed_finished=$(du */*final.txt | grep -P "^4\t" | wc -l)
failed=$(($failed_started + $failed_finished))
echo "Number of failed pipelines: $failed"

if [[ $failed != 0 ]]; then
  echo "Generate file containing failed pipelines: failed.txt"
  du */ | grep -P "^4\t" | cut -f 2 | sed 's/^.*\///g' > failed.txt
  du */*final.txt | grep -P "^4\t" | cut -f 2 | sed 's/\/.*$//g' >> failed.txt
else
  touch time.txt
  for i in slurm-*; do grep "Bundled pipelines ended" $i | sed 's/\[.*\] Bundled pipelines ended in //g' >> time.txt; done
  sort -k 1n -k 2n -o time.txt time.txt
  max_time=$(tail -n 1 time.txt)
  echo "Generate file containing times each bundled pipelines took: time.txt"
  echo "Longest time a pipeline took was: $max_time"
fi

if [[ $failed != 0 && $num_dir != 512 ]]; then
  echo "Generate file containing all pipelines that need to be rerun: rerun.txt"
  cat missing.txt failed.txt > rerun.txt
fi
