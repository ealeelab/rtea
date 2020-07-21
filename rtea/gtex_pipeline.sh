#!/bin/bash

<<DONE
# download hisat2 index
mkdir hg38
cd hg38
wget --content-disposition https://cloud.biohpc.swmed.edu/index.php/s/grch38_snp_tran/download
tar -xzf grch38_snp_tran.tar.gz
cd ..
DONE

# test run
ID="SRR811285"
R1="${ID}_1.fastq.gz"
R2="${ID}_2.fastq.gz"
HISAT_REF="hg38/grch38_snp_tran"
THREADS=8

./rtea_pipeline_without_scallop $R1 $R2 $ID $HISAT_REF $THREADS $PWD hg38 resume
