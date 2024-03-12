#!/bin/bash

# Define the source and destination directories
sourceDirectory="/path/to/sourceDirectory"
destinationDirectory="/path/to/destinationDirectory"

# Check if the source directory exists
if [[ ! -d "$sourceDirectory" ]]; then
    echo "Source directory '$sourceDirectory' does not exist."
    exit 1
fi

# Create the destination directory if it does not exist
if [[ ! -d "$destinationDirectory" ]]; then
    echo "Destination directory '$destinationDirectory' does not exists"
    exit 1
fi

# Loop through the subdirectories of sourceDirectory where subdirectories are sampleGroup
for sampleGroup in "$sourceDirectory"/*/; do
    # Get the name of the sampleGroup (basename of the directory path)
    sampleGroupName=$(basename "$sampleGroup")

    # Loop through subdirectories of sampleGroup which represent sampleID
    for sampleID in "$sampleGroup"*/; do
        # Get the name of the sampleID (basename of the directory path)
        sampleIDName=$(basename "$sampleID")

        # Create corresponding sampleGroup/sampleID directory in destinationDirectory if it doesn't exist
        targetDir="$destinationDirectory/$sampleGroupName/$sampleIDName"
        #mkdir -p "$targetDir"

        # Copy the *.rtea.depth.text files to corresponding directory in destinationDirectory
        filesToCopy=("$sampleID"*.rtea.depth.text)
        if [ ${#filesToCopy[@]} -gt 0 ]; then
            cp "${filesToCopy[@]}" "$targetDir"
        else
            echo "No files to copy for $sampleGroupName/$sampleIDName"
        fi
    done
done

echo "File copy process completed."
