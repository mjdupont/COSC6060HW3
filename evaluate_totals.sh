#!/bin/bash

mpicc -o DLBPolyEvaluation.exe DLBPolyEvaluation.c

rm -r -f "results"
mkdir "results"
for file in results/queue.csv; do
  echo "Procs,Terms,ChunkSize,Variable,Total,MaxTime" > $file
done

for polysize in 10 100 1000 10000 50000; do
  for np in 1 2 4 8 16; do
    for chunkSize in 1 5 10 20 50 100; do
      mpirun -np $np DLBPolyEvaluation.exe -t -p $polysize -c $chunkSize >> results/queue.csv
      echo "Completed $polysize terms with $np procs and $chunkSize chunkSize"
    done
  done
done

echo "Done."



