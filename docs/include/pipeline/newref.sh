#!/usr/bin/env bash

# PARAMETERS

CORES=8
INPUT_DIR="path/to/converts_ref" # folder containing reference .npz files
OUTPUT_DIR="path/to/references" # existing folder
REF_SIZES="1000 500 200 100" # space separated list (kb)

# SCRIPT

for REF in ${REF_SIZES}
do
    echo "Creating reference at bins size ${REF} kb"

    WisecondorX newref ${INPUT_DIR}/*.npz \
    ${OUTPUT_DIR}/reference.${REF}kb.npz \
    --binsize ${REF}000 --cpus ${CORES}
done