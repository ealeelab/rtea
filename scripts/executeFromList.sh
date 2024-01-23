#!/bin/bash

# Specify the file containing paths
file_path="/lab-share/Gene-Lee-e2/Public/share/projects/rTea/rtea-python/scripts/test.dirs"

# Check if the file exists
if [ ! -f "$file_path" ]; then
    echo "Error: File '$file_path' not found."
    exit 1
fi

# Use a while loop to iterate over each line in the file
while IFS= read -r line; do
    # Extract the components using parameter expansion
    full_line="$line"
    sampleName="${line##*/}"
    sampleGroup="${line%/*}"
    sampleGroup="${sampleGroup##*/}"


    # Print or use the variables as needed
    echo "Sample Location: $full_line"
    echo "Sample Group: $sampleGroup"
    echo "Sample Name: $sampleName"
done < "$file_path"
