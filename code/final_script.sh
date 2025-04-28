#!/bin/bash

# Define the list of files
FILES=(
  "RV11/BB11015.tfa"
  "RV11/BB11028.tfa"
  "RV11/BBS11015.tfa"
  "RV12/BB12004.tfa"
  "RV12/BBS12011.tfa"
)

LONGFEW=(
  "RV11/BB11016.tfa"
  "RV11/BB11026.tfa"
  "RV11/BB11037.tfa"
  "RV11/BB11038.tfa"
)
# Loop through each file
for file in "${LONGFEW[@]}"; do
    # Extract the filename (e.g., BB20001.tfa)
    base_filename=$(basename "$file" .tfa)

    # Paths
    input_path="../data/bb3_release/${file}"
    output_path="../data/output/${base_filename}.out"

    # Run sequential version and save output
    echo "Running bm_seq on ${file}"
    ./bm_seq -i "${input_path}" -o temp_seq_output.txt > temp_seq_terminal.txt 2>&1

    # Run parallel version and save output
    echo "Running bm_par on ${file}"
    mpirun -np 8 ./bm_par -i "${input_path}" -o temp_par_output.txt > temp_par_terminal.txt 2>&1

    # Combine the terminal outputs into a single .out file
    {
      echo "===== Sequential Run Output ====="
      cat temp_seq_terminal.txt
      echo ""
      echo "===== Parallel Run Output ====="
      cat temp_par_terminal.txt
    } > "${output_path}"

    # Clean up temporary files
    rm temp_seq_output.txt temp_par_output.txt temp_seq_terminal.txt temp_par_terminal.txt

    # Info
    echo "Processed $file -> ${output_path}"
done
