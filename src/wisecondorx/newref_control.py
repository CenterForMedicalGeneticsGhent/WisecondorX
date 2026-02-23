# WisecondorX

import logging
import os
import sys
import time
from typing import List, Optional

import numpy as np
from concurrent import futures

from wisecondorx.newref_tools import (
    normalize_and_mask,
    train_pca,
    get_reference,
)

"""
Outputs preparation files of read depth normalized
data and contains PCA information to execute between-
sample normalization during testing. Function is
executed three times. Once for autosomes, once for XX
gonosomes (if enough females are included) and once
for XY gonosomes (if enough males are included).
"""


def tool_newref_prep(
    prepdatafile: str,
    prepfile: str,
    binsize: int,
    samples: np.ndarray,
    gender: str,
    mask: np.ndarray,
    bins_per_chr: List[int],
) -> None:
    if gender == "A":
        last_chr = 22
    elif gender == "F":
        last_chr = 23
    else:
        last_chr = 24

    bins_per_chr = bins_per_chr[:last_chr]
    mask = mask[: np.sum(bins_per_chr)]

    masked_data = normalize_and_mask(samples, range(1, last_chr + 1), mask)
    pca_corrected_data, pca = train_pca(masked_data)

    masked_bins_per_chr = [
        sum(mask[sum(bins_per_chr[:i]) : sum(bins_per_chr[:i]) + x])
        for i, x in enumerate(bins_per_chr)
    ]
    masked_bins_per_chr_cum = [
        sum(masked_bins_per_chr[: x + 1])
        for x in range(len(masked_bins_per_chr))
    ]

    np.save(prepdatafile, pca_corrected_data)

    np.savez_compressed(
        prepfile,
        binsize=binsize,
        gender=gender,
        mask=mask,
        bins_per_chr=bins_per_chr,
        masked_bins_per_chr=masked_bins_per_chr,
        masked_bins_per_chr_cum=masked_bins_per_chr_cum,
        pca_components=pca.components_,
        pca_mean=pca.mean_,
    )


"""
Prepares subfiles if multi-threading is requested.
Main file is split in 'cpus' subfiles, each subfile
is processed by a separate thread.
"""


def tool_newref_main(
    prepdatafile: str,
    prepfile: str,
    partfile: str,
    tmpoutfile: str,
    refsize: int,
    cpus: int,
    part: Optional[List[int]] = None,
) -> None:
    pca_corrected_data = np.load(prepdatafile, mmap_mode="r")
    if cpus != 1:
        with futures.ThreadPoolExecutor(max_workers=cpus) as executor:
            for p in range(1, cpus + 1):
                executor.submit(
                    _tool_newref_part,
                    prepfile,
                    partfile,
                    refsize,
                    [p, cpus],
                    pca_corrected_data,
                )
            executor.shutdown(wait=True)
    else:
        for p in range(1, cpus + 1):
            _tool_newref_part(
                prepfile, partfile, refsize, [p, cpus], pca_corrected_data
            )

    tool_newref_post(prepfile, partfile, tmpoutfile, cpus)

    os.remove(prepfile)
    os.remove(prepdatafile)
    for p in range(1, cpus + 1):
        os.remove("{}_{}.npz".format(partfile, str(p)))


"""
Function executed once for each thread. Controls
within-sample reference creation.
"""


def _tool_newref_part(
    prepfile: str,
    partfile: str,
    refsize: int,
    part: List[int],
    pca_corrected_data: np.ndarray,
) -> None:
    if part[0] > part[1]:
        logging.critical(
            "Part should be smaller or equal to total parts:{} > {} is wrong".format(
                part[0], part[1]
            )
        )
        sys.exit()
    if part[0] < 0:
        logging.critical(
            "Part should be at least zero: {} < 0 is wrong".format(part[0])
        )
        sys.exit()

    npzdata = np.load(prepfile, encoding="latin1", allow_pickle=True)
    masked_bins_per_chr = npzdata["masked_bins_per_chr"]
    masked_bins_per_chr_cum = npzdata["masked_bins_per_chr_cum"]

    indexes, distances, null_ratios = get_reference(
        pca_corrected_data,
        masked_bins_per_chr,
        masked_bins_per_chr_cum,
        ref_size=refsize,
        part=part[0],
        split_parts=part[1],
    )

    np.savez_compressed(
        "{}_{}.npz".format(partfile, str(part[0])),
        indexes=indexes,
        distances=distances,
        null_ratios=null_ratios,
    )


"""
Merges separate subfiles (one for each thread) to a
new temporary output file.
"""


def tool_newref_post(
    prepfile: str, partfile: str, tmpoutfile: str, cpus: int
) -> None:
    npzdata_prep = np.load(prepfile, encoding="latin1", allow_pickle=True)

    big_indexes = []
    big_distances = []
    big_null_ratios = []
    for p in range(1, cpus + 1):
        infile = "{}_{}.npz".format(partfile, str(p))
        npzdata_part = np.load(infile, encoding="latin1")
        big_indexes.extend(npzdata_part["indexes"])
        big_distances.extend(npzdata_part["distances"])
        big_null_ratios.extend(npzdata_part["null_ratios"])

    indexes = np.array(big_indexes)
    distances = np.array(big_distances)
    null_ratios = np.array(big_null_ratios)

    np.savez_compressed(
        tmpoutfile,
        binsize=npzdata_prep["binsize"].item(),
        gender=npzdata_prep["gender"].item(),
        mask=npzdata_prep["mask"],
        bins_per_chr=npzdata_prep["bins_per_chr"],
        masked_bins_per_chr=npzdata_prep["masked_bins_per_chr"],
        masked_bins_per_chr_cum=npzdata_prep["masked_bins_per_chr_cum"],
        pca_components=npzdata_prep["pca_components"],
        pca_mean=npzdata_prep["pca_mean"],
        indexes=indexes,
        distances=distances,
        null_ratios=null_ratios,
    )


"""
Tries to remove text file, when it is busy, until becomes successful.
This function, prevents OSError: [Errno 26] Text file busy...
"""


def force_remove(file_id: str) -> None:
    attemp = 1
    while True:
        try:
            os.remove(file_id)
            break
        except Exception:
            print(
                "Attemp #{}: Cannot remove {}, because it is busy, trying again...".format(
                    attemp, file_id
                )
            )
            attemp = attemp + 1
            time.sleep(5)


"""
Merges separate subfiles (A, F, M) to one final
reference file.
"""


def tool_newref_merge(
    outfile: str, nipt: bool, outfiles: List[str], trained_cutoff: float
) -> None:
    final_ref = {"has_female": False, "has_male": False}
    for file_id in outfiles:
        npz_file = np.load(file_id, encoding="latin1", allow_pickle=True)
        gender = str(npz_file["gender"])
        for component in [x for x in npz_file.keys() if x != "gender"]:
            if gender == "F":
                final_ref["has_female"] = True
                final_ref["{}.F".format(str(component))] = npz_file[component]
            elif gender == "M":
                final_ref["has_male"] = True
                final_ref["{}.M".format(str(component))] = npz_file[component]
            else:
                final_ref[str(component)] = npz_file[component]
        force_remove(file_id)
    final_ref["is_nipt"] = nipt
    final_ref["trained_cutoff"] = trained_cutoff
    np.savez_compressed(outfile, **final_ref)
