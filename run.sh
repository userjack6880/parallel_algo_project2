#!/bin/bash
for run in 1 2 10 25 50 100
do
  # no increase
  for threads in 01 02 04 08 16 32 64
  do
    sbatch --output=/dev/null --job-name="run $run $threads threads" run$threads.js $run
  done
done