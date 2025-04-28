#!/bin/bash

# for sequential run
./bm_seq -i "../data/bb3_release/RV20/BB20001.tfa" -o out_seq.txt


# for parallel run

mpirun -np 8 ./bm_par -i "../data/bb3_release/RV20/BB20001.tfa" -o out_par.txt

# Script shell for running 