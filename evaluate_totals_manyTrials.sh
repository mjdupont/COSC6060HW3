#!/bin/bash

mpicc -o DLBPolyEvaluation.exe DLBPolyEvaluation.c

rm -r -f "results"
mkdir "results"
for file in results/queue.csv; do
  echo "Procs,Terms,ChunkSize,Variable,Total,MaxTime" > $file
done

for trial in $(seq 1 16); do
  for polysize in 50000 75000; do
    for np in 2 4 8 16; do
      for chunkSize in 1 2 4 8 16 32; do
        mpirun -np $np DLBPolyEvaluation.exe --terse -p $polysize -c $chunkSize >> results/queue.csv
      done
    done
  done
done

echo "Done."


