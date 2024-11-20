#!/bin/bash

# Bash script to get median and max SVlen from vcf file
# Check if VCF file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <vcf_file>"
    exit 1
fi

vcf_file=$1

# Extract SVLEN values, convert negative values to positive (absolute value), sort, and calculate median and max
svlen_values=$(grep -v "^#" "$vcf_file" | \
    awk '{
        # Loop through all INFO fields and look for SVLEN=
        split($8, info, ";")
        for (i in info) {
            if (info[i] ~ /^SVLEN=/) {
                # Extract the SVLEN value, remove "SVLEN=", and convert to absolute value
                svlen_value = substr(info[i], 7)
                if (svlen_value < 0) {
                    svlen_value = -svlen_value
                }
                print svlen_value
            }
        }
    }' | sort -n)

# Calculate median
median_svlen=$(echo "$svlen_values" | awk '{arr[NR]=$1} END {if (NR % 2) {print arr[(NR+1)/2]} else {print (arr[NR/2] + arr[NR/2+1])/2}}')

# Calculate max
max_svlen=$(echo "$svlen_values" | tail -n 1)

# Output the results
echo "Median SVLEN: $median_svlen"
echo "Max SVLEN: $max_svlen"
