#!/bin/bash

# Specify the file containing paths
# Start Conda Env

$CONDA_HOME=/home/aleelab/miniforge3
. $CONDA_HOME/etc/profile.d/conda.sh

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file_path>"
    exit 1
fi

file_path="$1"

# Check if the destination directory exists
if [ ! -d "$file_path" ]; then
    echo "Error: Directory '$file_path' not found."
    exit 1
fi

#file_path="/lab-share/Gene-Lee-e2/Public/share/projects/rTea/rtea-python/scripts/test.dirs"

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

    # Program Start

done < "$file_path"
