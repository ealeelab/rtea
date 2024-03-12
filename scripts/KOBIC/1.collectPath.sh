#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <destination_directory>"
    exit 1
fi

destination_directory="$1"

# Check if the destination directory exists
if [ ! -d "$destination_directory" ]; then
    echo "Error: Directory '$destination_directory' not found."
    exit 1
fi

# Get real path of the destination directory
destination_directory=$(realpath "$destination_directory")

# Collect real paths of directories at depth 2
directories_at_depth_2=$(find "$destination_directory" -maxdepth 2 -mindepth 2 -type d -exec realpath {} \;)

# Print the result
#echo "Real paths of directories at depth 2:"
echo "$directories_at_depth_2"

