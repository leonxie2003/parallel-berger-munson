#!/bin/bash

INPUT_FILES=(
    "demo/few_long.tfa"
    "demo/few_very_long.tfa"
)

mkdir "../data/output"

for input_file in "${INPUT_FILES[@]}"; do
    echo "=========="
    base_filename=$(basename "$input_file" .tfa)

    input_path="../data/${input_file}"
    output_path="../data/output/${base_filename}"

    echo "Running bm_seq on ${input_file}"
    ./bm_seq -i "${input_path}" -o "${output_path}_seq.out" -r P
    echo ""

    echo "Running bm_par (p=2) on ${input_file}"
    mpirun -np 2 ./bm_par -i "${input_path}" -o "${output_path}_par_2.out" -r P
    echo ""

    echo "Running bm_par (p=4) on ${input_file}"
    mpirun -np 4 ./bm_par -i "${input_path}" -o "${output_path}_par_4.out" -r P
    echo ""

    echo "Running bm_par (p=8) on ${input_file}"
    mpirun -np 8 ./bm_par -i "${input_path}" -o "${output_path}_par_8.out" -r P

    echo "=========="
done
