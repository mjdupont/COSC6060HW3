#!/bin/bash

mpicc -o DLBPolyEvaluation.exe DLBPolyEvaluation.c

rm -r -f "results_individual"
mkdir "results_individual"
for file in results_individual/queue.csv; do
  echo "Procs,Terms,ChunkSize,Variable,Total,Rank,Time" > $file
done

for polysize in 50000 75000; do
  for np in 2 4 8 16; do
    for chunkSize in 1 2 4 8 16 32; do
      mpirun -np $np DLBPolyEvaluation.exe --veryTerse -p $polysize -c $chunkSize >> results_individual/queue.csv
    done
  done
done

echo "Done."