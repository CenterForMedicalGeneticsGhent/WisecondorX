# This script reformats samples to make them backwards compatible between WisecondorX versions.
# It can also be used to transform .npz files from the original WISECONDOR to WisecondorX.
# Note that is only compatible with .npz files resulting from the convert version, other
# functions (newref, predict) should be re-run.

# Usage: python fix.convert.npz.py input.npz output.npz

import sys
import numpy as np

_, in_npz, out_npz = sys.argv

def toolConvert():
    1
    
def get_gender(sample):

    tot_reads = float(sum([sum(sample[str(x)]) for x in range(1, 25)]))
    X_reads = float(sum(sample["23"]))
    X_len = float(len(sample["23"]))
    Y_reads = float(sum(sample["24"]))
    Y_len = float(len(sample["24"]))

    X = (X_reads / tot_reads) / X_len * 0.5
    Y = (Y_reads / tot_reads) / Y_len

    # X/Y               = ?
    # 1/1 (MALE)        = 1
    # 2/noise (FEMALE)  = [4,8]
    # cut-off 3 -- should be robust vs noise and mosaic large subchromosomal duplication/deletions
    if X/Y < 3:
        return "M"
    else:
        return "F"

def reformat(sample):
    for chr in sample.keys():
        data = sample[chr]
        if chr == "X":
            chr = "23"
            sample.pop("X")
        if chr == "Y":
            chr = "24"
            sample.pop("Y")
        sample[chr] = data
    return sample

npz = np.load(in_npz)
sample = npz["sample"].item()
sample = reformat(sample)
gender = get_gender(sample)

if "binsize" in npz.keys():
    binsize = npz["binsize"].item()
else:
    binsize = npz["arguments"].item()["binsize"]

np.savez_compressed(out_npz,
    binsize=binsize,
    sample=sample,
    gender=gender,
    quality=npz["quality"].item())

print("Succes!")
