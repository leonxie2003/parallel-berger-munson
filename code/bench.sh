#!/bin/bash

INPUT_FILES=(
    # few short
    # "RV11/BB11013.tfa" # 5 x 115

    # few long
    # "RV12/BB12028.tfa" # 8 x 519
    # "RV12/BB12007.tfa" # 8 x 1263
    # "RV30/BB30027.tfa" # 20 x 239

    # many short
    # "RV30/BBS30009.tfa" # 38 x 74
    # "RV20/BBS20030.tfa" # 47 x 97
    # "RV20/BBS20032.tfa" # 61 x 116

    # one that didn't scale well
    "RV20/BB20007.tfa" # 23 x 515
)

for input_file in "${INPUT_FILES[@]}"; do
    echo "=========="
    base_filename=$(basename "$input_file" .tfa)

    input_path="../data/bb3_release/${input_file}"
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
