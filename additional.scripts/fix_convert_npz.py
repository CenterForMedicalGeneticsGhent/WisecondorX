#!/usr/bin/env python

import argparse
import numpy as np

def toolConvert():
    # Added for compatibility with with original WiSeCondor files
    return True

def get_gender(sample):
    tot_reads = float(sum([sum(sample[str(x)]) for x in range(1, 25)]))
    x_reads = float(sum(sample["23"]))
    x_len = float(len(sample["23"]))
    y_reads = float(sum(sample["24"]))
    y_len = float(len(sample["24"]))

    X = (x_reads / tot_reads) / x_len * 0.5
    Y = (y_reads / tot_reads) / y_len

    # X/Y               = ?
    # 1/1 (MALE)        = 1
    # 2/noise (FEMALE)  = [4,8]
    # cut-off 3 -- should be robust vs noise and mosaic large subchromosomal duplication/deletions
    if X/Y < 3:
        return "M"
    else:
        return "F"


def reformat(sample):
    for chrom in sample.keys():
        data = sample[chrom]
        if chrom == "X":
            chrom = "23"
            sample.pop("X")
        if chrom == "Y":
            chrom = "24"
            sample.pop("Y")
        sample[chrom] = data
    return sample


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="\
            This script reformats samples to make them backwards compatible between WisecondorX versions. \
            It can also be used to transform .npz files from the original WISECONDOR to WisecondorX. \
            Note that is only compatible with .npz files resulting from the convert version, \
            other functions (newref, predict) should be re-run."
                                     )
    parser.add_argument('in_npz', help='Input npz from Wisecondor or WisecondorX <0.2.0')
    parser.add_argument('out_npz', help='Output npz')

    args = parser.parse_args()

    npz = np.load(args.in_npz)
    sample = npz["sample"].item()
    sample = reformat(sample)
    gender = get_gender(sample)

    if "binsize" in npz.keys():
        binsize = npz["binsize"].item()
    else:
        binsize = npz["arguments"].item()["binsize"]

    np.savez_compressed(args.out_npz, binsize=binsize, sample=sample, gender=gender, quality=npz["quality"].item())
    print("Succes!")
