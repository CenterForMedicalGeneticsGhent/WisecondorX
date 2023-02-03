#!/usr/bin/env bash

# PARAMETERS

NPZ_FILES="path/to/npz_files.txt" # file
# Example of the structure of this file:
# ID_1 path/to/converts_test/ID_1.npz
# ID_2 path/to/converts_test/ID_2.npz
# ...
REF="path/to/references/reference.100kb.npz"
OUTPUT_DIR="path/to/output_100kb" # existing folder

# SCRIPT

while read LINE; do

    SAMPLE=$(echo $LINE | awk -F ' ' '{print $1}')
    NPZ=$(echo $LINE | awk -F ' ' '{print $2}')

    echo "Predicting sample ${SAMPLE}"
    WisecondorX predict ${NPZ} ${REF} ${OUTPUT_DIR}/${SAMPLE} --plot --bed

done <${NPZ_FILES}