#!/bin/bash

# Specify the directory containing text files
directory="/path/to/your/directory"
programLocation=/home/aleelab/rtea

# Check if the directory exists
if [ -d "$directory" ]; then
    # Iterate over each text file in the directory
    for file in "$directory"/*.txt; do
        # Check if the file is a regular file
        if [ -f "$file" ]; then
            # Perform actions on the text file, for example, print its content
            echo "Processing file: $file"
            cat "$file"  # Replace this with your actual processing logic

            $programLocation/scripts/executeFromList.sh $file

        fi
    done
else
    echo "Directory does not exist: $directory"
fi
