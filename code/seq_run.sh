echo "Hi"

# List of files you want to process
FILES=(
 # "RV12/BBS12014.tfa
 # "RV40/BB40007.tfa"
 # "RV30/BB30029.tfa"
  "RV20/BBS20036.tfa"
 # "RV30/BB30003.tfa"
)

# Loop through each file
for file in "${FILES[@]}"; do
    # Get just the filename without folder (basename)
    filename="$file"


    # Run the command
    ./bm_seq -i ../data/bb3_release/${filename} -o ${filename}.txt

    # Print info
    echo "Processed $file -> ${filename}.txt"
done
