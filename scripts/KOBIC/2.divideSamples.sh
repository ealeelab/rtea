#!/bin/bash
# Specify the input file

# Check if the correct number of arguments is provided
output_file_prefix=tcga_

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sample list>"
    exit 1
fi

input_file="$1"


echo $input_file

# Check if the file exists
if [ ! -f $input_file ]; then
    echo “Error: File ‘$input_file’ not found.”
    exit 1
fi

# Calculate the number of lines per file
total_lines=$(wc -l < $input_file)
lines_per_file=$((total_lines / 20))
remainder=$((total_lines % 20))

echo "total line: $total_lines"
echo "lines_pser_file: $lines_per_file"
echo "remainder: $remainder"

# Use the split command to divide the file
if [ $lines_per_file -gt 0 ]; then
    split -l "$lines_per_file" "$input_file" $output_file_prefix
else
    echo "Error: The file is empty or has fewer than 20 lines."
    exit 1
fi

# Handle the remainder lines
if [ $remainder -gt 0 ]; then
    tail -n "$remainder" "$input_file" >> "${output_file_prefix}aa"
fi

echo "File divided into 20 parts with remainder lines included."