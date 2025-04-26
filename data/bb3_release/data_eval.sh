
# Set the base folder (first argument or default to ".")
BASE_FOLDER="${1:-.}"

# Output file
OUTPUT_FILE="msf_analysis_report.txt"

# Empty the output file
> "$OUTPUT_FILE"

# Find all .msf files
find "$BASE_FOLDER" -type f -name "*.msf" | while read -r file; do
    {
    echo "Processing: $file"

    # Extract and process "Name:" lines with "Len:"
    awk '
    BEGIN { total_len=0; count=0 }
    /^ Name:/ {
        if ($0 ~ /Len:[ ]*[0-9]+/) {
            match($0, /Len:[ ]*([0-9]+)/, arr)
            if (arr[1] != "") {
                total_len += arr[1]
                count++
            }
        }
    }
    END {
        if (count > 0) {
            avg_len = total_len / count
            printf "Number of sequences: %d\n", count
            printf "Average sequence length: %.2f\n", avg_len
        } else {
            print "No sequences found!"
        }
    }
    ' "$file"

    echo "-------------------------------"
    } >> "$OUTPUT_FILE"
done

echo "Report saved to $OUTPUT_FILE"
