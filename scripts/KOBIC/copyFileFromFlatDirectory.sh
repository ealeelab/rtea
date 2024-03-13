#!/bin/bash

# Main directory containing subdirectories with sample files
#main_dir="/BiO/data/PROJECT/HMS/tcga/ESCA-recovery"
main_dir="/BiO/data/PROJECT/HMS/TCGA/ESCA"
data_dir="/BiO/data/PROJECT/HMS/TCGA/ESCA-recovery"

# Move into the main directory
cd "$data_dir" || exit

# Recursively find files with sample IDs in subdirectories
find $main_dir -type d -mindepth 1 -maxdepth 1 -print0 | while IFS= read -r -d '' directory; do
    # Extract sample ID from file name
    sample_id="${directory##*/}"

    echo $sample_id
    echo $directory


    # Loop over files starting with "sampleID"
    for file in ${sample_id}*; do
        # Extract sample ID from file name
        # sample_id=$(echo "$file" | cut -d '_' -f 1)

        if [[ $file == ${sample_id} ]]; then
            continue
        fi

        if [[ -f $file ]]; then
           echo $file
           mv "$file" "$directory/"
        fi

        #echo ""
        # Create directory if it doesn't exist
        #mkdir -p "$sample_id"

        # Move file into directory
        #mv "$file" "$directory/"

    done
done

