#!/usr/bin/env bash

# PARAMETERS

BAM_FILES="path/to/bam_files.txt" # file
# Example of the structure of this file:
# ID_1 path/to/ID_1.bam
# ID_2 path/to/ID_2.bam
# ...
OUTPUT_DIR="path/to/converts" # existing folder

# SCRIPT

while read LINE; do

    SAMPLE=$(echo $LINE | awk -F ' ' '{print $1}')
    BAM=$(echo $LINE | awk -F ' ' '{print $2}')

    echo "Creating 5kb bins for sample ${SAMPLE}"
    WisecondorX convert ${BAM} ${OUTPUT_DIR}/${SAMPLE}.npz

done < ${BAM_FILES}