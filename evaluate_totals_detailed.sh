#!/bin/bash

mpicc -o DLBPolyEvaluation.exe DLBPolyEvaluation.c

rm -r -f "results_detailed"
mkdir "results_detailed"
for file in results/queue.csv; do
  echo "Procs,Terms,ChunkSize,Variable,Total,MaxTime" > $file
done

for polysize in 1000; do
  for np in 1 2 3 5 9 13 16; do
    for chunkSize in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 50 100; do
      mpirun -np $np DLBPolyEvaluation.exe -t -p $polysize -c $chunkSize >> results/queue.csv
      echo "Completed $polysize terms with $np procs and $chunkSize chunkSize"
    done
  done
done

echo "Done."



