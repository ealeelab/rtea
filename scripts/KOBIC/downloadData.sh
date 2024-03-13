#!/bin/bash

#$ -q kobic.q
#$ -N download
#$ -e /home/aleelab/logs
#$ -o /home/aleelab/logs
#$ -S /bin/bash
#$ -cwd


CELLZENTEK=220.125.199.66

#DATA_PATH=/data/data2_old/boramlee/tei/rtea/bam
DATA_PATH=/data/data2_old/boramlee/tei/rtea/rtea

set -x

rsync -avz -e "ssh -i ~/.ssh/id_cellzentek" --progress aleelab@${CELLZENTEK}:${DATA_PATH} /BiO/scratch/users/aleelab/copm_rtea
