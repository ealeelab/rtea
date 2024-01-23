#!/bin/bash
# Specify the input file
input_file=“/path/to/your/input_file.txt”
# Check if the file exists
if [ ! -f “$input_file” ]; then
    echo “Error: File ‘$input_file’ not found.”
    exit 1
fi

# Calculate the number of lines per file
total_lines=$(wc -l < “$input_file”)
lines_per_file=$((total_lines / 20))
remainder=$((total_lines % 20))

# Use the split command to divide the file
split -l “$lines_per_file” “$input_file” “output_file_prefix”
# Handle the remainder lines
if [ “$remainder” -gt 0 ]; then
    tail -n “$remainder” “$input_file” >> “output_file_prefixaa”
fi
echo “File divided into 20 parts with remainder lines included.”