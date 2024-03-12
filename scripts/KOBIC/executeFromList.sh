#!/bin/bash
#$ -q kobic.q
#$ -N copm_depth
#$ -e /BiO/scratch/users/aleelab/logs
#$ -o /BiO/scratch/users/aleelab/logs
#$ -S /bin/bash
#$ -pe pe_slots 4
#$ -R yes
#$ -l h_rt=23:59:59
#$ -l mem_free=16G
#$ -l h_vmem=60G
#$ -M junseokpark.kr@gmail.com
#$ -m e
#$ -cwd

# Specify the file containing paths
# Start Conda Env

CONDA_HOME=/home/aleelab/miniforge3
. $CONDA_HOME/etc/profile.d/conda.sh

conda activate rtea

PATH=/home/aleelab/rtea/bin:$PATH
SAVED_PATH=/BiO/scratch/users/aleelab

REPO_NAME="COPM"
CORES=4

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file_path>"
    exit 1
fi

file_path="$1"

# Check if the destination directory exists
if [ ! -f "$file_path" ]; then
    echo "Error: Directory '$file_path' not found."
    exit 1
fi

#file_path="/lab-share/Gene-Lee-e2/Public/share/projects/rTea/rtea-python/scripts/test.dirs"

# Use a while loop to iterate over each line in the file
while IFS= read -r line; do
    # Extract the components using parameter expansion
    full_line="$line"
    sampleName="${line##*/}"
    sampleName=$(echo "$sampleName" | cut -d'.' -f1)


    sampleGroup="${line%/*}"
    sampleGroup="${sampleGroup##*/}"

    # Print or use the variables as needed
    echo "Sample Location: $full_line"
    echo "Sample Group: $sampleGroup"
    echo "Sample Name: $sampleName"

    # Create directory
    OUTPUT_DIR="${SAVED_PATH}/${REPO_NAME}/${sampleGroup}/${sampleName}"
    OUTPUT_DIR="${SAVED_PATH}/${REPO_NAME}/${sampleGroup}"
    echo "Output Dir: $OUTPUT_DIR"

    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Create output directory"
        mkdir -p $OUTPUT_DIR
    else
        echo "there is exising directory. It will be removed then we will re-create it"
        #rm -rf $OUTPUT_DIR
        #mkdir -p $OUTPUT_DIR
    fi 

    # Program Start
    line=$(dirname "$line")

    if [ -f "${OUTPUT_DIR}/${sampleName}.hisat2.rtea.depth.txt" ]; then
        echo "The file ${sampleName}.hisat2.rtea.depth.txt exists. Continuing with further actions."
        continue
    else
        echo "The file ${sampleName}.hisat2.rtea.depth.txt does not exist. Exiting script."
        depthCalculation -r ${line}/${sampleName}.rtea.txt -b ${line}/${sampleName}.hisat2.bam -o ${OUTPUT_DIR}/${sampleName}.hisat2.rtea.depth.txt -c $CORES
    fi

done < "$file_path"
