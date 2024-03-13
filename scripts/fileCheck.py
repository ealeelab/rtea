import os

def count_txt_files(directory):
    txt_file_count = 0
    bam_file_count = 0
    # Recursively walk through all directories and subdirectories
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".txt"):
                txt_file_count += 1
            if file.endswith('.bam'):
                bam_file_count += 1

    return txt_file_count, bam_file_count

# Specify the parent directory
parent_directory = '/path/to/parent/directory'  # Replace with the path to your parent directory
parent_directory = '/BiO/data/PROJECT/HMS/temp/icgc'
#parent_directory = '/BiO/data/PROJECT/HMS/tcga'

# Count the number of .txt files
num_txt_files, num_bam_files = count_txt_files(parent_directory)

# Print the count
print("Number of .txt files:", num_txt_files)
print("Number of .bam files:", num_bam_files)
