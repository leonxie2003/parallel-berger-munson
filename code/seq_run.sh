echo "Hi"

# List of files you want to process
# FILES=(
#  # "RV12/BBS12014.tfa
#  # "RV40/BB40007.tfa"
#  # "RV30/BB30029.tfa"
#   "RV20/BBS20036.tfa"
#  # "RV30/BB30003.tfa"
# )

# FILES2=(
#   "RV11/BB11001.tfa"
#   "RV11/BBS11002.tfa"
#   "RV12/BB12002.tfa"
#   "RV12/BB12004.tfa"
#   "RV12/BB12008.tfa"
# )

# Around similar as BF12008, where I expect around 100 seconds to run
FILES=(
  "RV12/BBS12014.tfa"
  "RV20/BBS20036.tfa"
  "RV20/BB20001.tfa"
  "RV12/BB12013.tfa"
  "RV12/BBS12011.tfa"
  "RV20/BB20007.tfa"
)


# # Similar as BF12008, around twice longer time expected
# FILES4=(
#   "RV20/BB20004.tfa"
#   "RV20/BB20006.tfa"
#   "RV20/BBS20006.tfa"
# )

# Loop through each file
for file in "${FILES[@]}"; do
    # Get just the filename without folder (basename)
    filename="$file"


    # Run the command
    ./bm_seq -i ../data/bb3_release/${filename} -o ${filename}.txt

    # Print info
    echo "Processed $file -> ${filename}.txt"
done
