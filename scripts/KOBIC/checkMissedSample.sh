#!/bin/bash

# Define the source and destination directories
sourceDirectory="/BiO/data/PROJECT/HMS/tcga"

# Check if the source directory exists
if [[ ! -d "$sourceDirectory" ]]; then
    echo "Source directory '$sourceDirectory' does not exist."
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

        # Copy the *.rtea.depth.text files to corresponding directory in destinationDirectory
        fileToCheck=${sourceDirectory}/${sampleGroupName}/${sampleIDName}/${sampleIDName}.hisat2.rtea.depth.txt
        if [ -f $fileToCheck ]; then
            if [ -f ${sourceDirectory}/${sampleGroupName}/${sampleIDName}/${sampleIDName}.hisat2.bam ]; then
                echo "${sourceDirectory}/${sampleGroupName}/${sampleIDName}/${sampleIDName}.hisat2.bam is not exist"
            else
                echo ${sourceDirectory}/${sampleGroupName}/${sampleIDName}
            fi
        fi
    done
done