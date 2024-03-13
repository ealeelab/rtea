#!/bin/bash

# Main directory containing subdirectories with sample files
#main_dir="/BiO/data/PROJECT/HMS/tcga/ESCA-recovery"
main_dir="/BiO/data/PROJECT/HMS/icgc/res_rtea"

# Move into the main directory
cd "$main_dir" || exit

# Recursively find files with sample IDs in subdirectories
find . -type f -mindepth 1 -maxdepth 1 -name '*.bam' -print0 | while IFS= read -r -d '' file; do
    # Extract sample ID from file name
    sample_id=$(basename "$file" | cut -d '.' -f 1)

    echo $sample_id


    # Loop over files starting with "sampleID"
    for file in ${sample_id}*; do
        # Extract sample ID from file name
        # sample_id=$(echo "$file" | cut -d '_' -f 1)

        if [[ $file == ${sample_id} ]]; then
            continue
        fi

        echo $file
        #echo ""
        # Create directory if it doesn't exist
        mkdir -p "$sample_id"

        # Move file into directory
        mv "$file" "$sample_id/"

    done
done
