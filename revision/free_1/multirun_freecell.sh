#!/bin/bash

declare -a parameters=(
  0.03
  0.1
  0.3
  0.8
  1.0
  1.3
  1.7
  2.0
  3.0
  5.0
  10.0
)

for i in "${parameters[@]}"
do
  echo "[SHELL] Running test with parameter $i";
  nvcc Hazard.cu -o free_cell_${i}.out && ./free_cell_${i}.out $i;
  echo "[SHELL] test with parameter $i is done."
done
