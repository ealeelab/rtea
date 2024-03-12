#!/bin/bash


# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <located_directory> <program_location> <prefix>"
    exit 1
fi

directory="$1"
#programLocation=/home/aleelab/rtea
programLocation="$2"
prefix="$3"


# Check if the destination directory exists
if [ ! -d $directory ]; then
    echo "Error: Directory '$directory' not found."
    exit 1
fi

# Check if the directory exists
if [ -d $directory ]; then
    # Iterate over each text file in the directory
    echo "sample directory: $directory"
    for file in $directory/${prefix}_*; do
        # Check if the file is a regular file
        if [ -f "$file" ]; then
            # Perform actions on the text file, for example, print its content
            echo "Processing file: $file"
            #cat "$file"  # Replace this with your actual processing logic

            qsub $programLocation/scripts/KOBIC/executeFromList.sh $file

        fi
    done
else
    echo "Directory does not exist: $directory"
fi

