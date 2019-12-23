#!/bin/bash

export LD_LIBRARY_PATH=/usr/local/bin/lib/lib:/usr/local/lib
export ctea_base=/usr/local/bin/lib/ctea

BAM=$1
OUT_DIR=$2
BAM_FILE=${OUT_DIR}/$(basename $BAM)
CTEA_EXEC="/usr/local/bin/ctea_exec"
BWA_FA_FILE="/usr/local/bin/res/repeat_LINE1_ALU_SVA_HERV_human_youngTE.fa"

rm -f ${BAM_FILE}
rm -f ${BAM_FILE}.bai
ln -s $(realpath $BAM) ${BAM_FILE}
ln -s $(realpath ${BAM}.bai) ${BAM_FILE}.bai

${CTEA_EXEC} cteax -b ${BAM_FILE} -r ${BWA_FA_FILE} --nofilter

