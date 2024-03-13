#!/bin/bash

# Define the source and destination directories
sourceDirectory="/BiO/data/PROJECT/HMS/ICGC"

# Check if the source directory exists
if [[ ! -d "$sourceDirectory" ]]; then
    echo "Source directory '$sourceDirectory' does not exist."
    exit 1
fi

# Loop through the subdirectories of sourceDirectory where subdirectories are sampleGroup
#for sampleGroup in "$sourceDirectory"/*/; do
    # Get the name of the sampleGroup (basename of the directory path)
#    sampleGroupName=$(basename "$sampleGroup")

    # Loop through subdirectories of sampleGroup which represent sampleID
    #for sampleID in "$sampleGroup"*/; do
    for sampleID in "$sourceDirectory"/*/; do #COPM and ICGC
        # Get the name of the sampleID (basename of the directory path)
        sampleIDName=$(basename "$sampleID")

        # Create corresponding sampleGroup/sampleID directory in destinationDirectory if it doesn't exist
        targetDir="$sourceDirectory/$sampleGroupName/$sampleIDName"
        
        # Full path of bam file
        bam_file="$targetDir/${sampleIDName}.hisat2.bam"

        # Full paths of rtea.txt and rtea.depth.txt files
        rtea_txt="$targetDir/${sampleIDName}.rtea.txt"
        rtea_depth_txt="$targetDir/${sampleIDName}.hisat2.rtea.depth.txt"

            # Check existence of bam fil
		if [ -f "$bam_file" ]; then
			bam_existence=1
		else
			bam_existence=0
		fi

		# Check existence of rtea.txt file
		if [ -f "$rtea_txt" ]; then
		    rtea_txt_existence=1
		else
		    rtea_txt_existence=0
		fi

		# Check existence of rtea.depth.txt file
		if [ -f "$rtea_depth_txt" ]; then
		    rtea_depth_txt_existence=1
		else
		    rtea_depth_txt_existence=0
        fi
        echo -e "$targetDir\t$sampleGroupName\t$sampleIDName\t$bam_existence\t$rtea_txt_existence\t$rtea_depth_txt_existence"
    done
#done

