echo "Hi"

# List of files you want to process
FILES=(
  "RV30/BB30003.msf"
  "RV20/BBS20036.msf"
  "RV30/BB30029.msf"
  "RV12/BBS12014.msf"
  "RV40/BB40007.msf"
)

# Loop through each file
for file in "${FILES[@]}"; do
    # Get just the filename without folder (basename)
    filename=$(basename "$file" .msf)
    
    # Run the command
    ./bm_seq -i "$file" -o "${filename}.txt"
    
    # Print info
    echo "Processed $file -> ${filename}.txt"
done