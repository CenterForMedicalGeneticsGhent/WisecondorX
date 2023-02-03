# WisecondorX

import copy
import logging
import os
import sys
import time

import numpy as np
from concurrent import futures

from wisecondorX.newref_tools import normalize_and_mask, train_pca, get_reference

"""
Outputs preparation files of read depth normalized
data and contains PCA information to execute between-
sample normalization during testing. Function is
executed three times. Once for autosomes, once for XX
gonosomes (if enough females are included) and once
for XY gonosomes (if enough males are included).
"""


def tool_newref_prep(args, samples, gender, mask, bins_per_chr):
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
        sum(masked_bins_per_chr[: x + 1]) for x in range(len(masked_bins_per_chr))
    ]

    np.savez_compressed(
        args.prepfile,
        binsize=args.binsize,
        gender=gender,
        mask=mask,
        masked_data=masked_data,
        bins_per_chr=bins_per_chr,
        masked_bins_per_chr=masked_bins_per_chr,
        masked_bins_per_chr_cum=masked_bins_per_chr_cum,
        pca_corrected_data=pca_corrected_data,
        pca_components=pca.components_,
        pca_mean=pca.mean_,
    )


"""
Prepares subfiles if multi-threading is requested.
Main file is split in 'cpus' subfiles, each subfile
is processed by a separate thread.
"""


def tool_newref_main(args, cpus):
    if cpus != 1:
        with futures.ThreadPoolExecutor(max_workers=args.cpus) as executor:
            for part in range(1, cpus + 1):
                this_args = copy.copy(args)
                this_args.part = [part, cpus]
                executor.submit(_tool_newref_part, this_args)
            executor.shutdown(wait=True)
    else:
        for part in range(1, cpus + 1):
            args.part = [part, cpus]
            _tool_newref_part(args)

    tool_newref_post(args, cpus)

    os.remove(args.prepfile)
    for part in range(1, cpus + 1):
        os.remove("{}_{}.npz".format(args.partfile, str(part)))


"""
Function executed once for each thread. Controls
within-sample reference creation.
"""


def _tool_newref_part(args):
    if args.part[0] > args.part[1]:
        logging.critical(
            "Part should be smaller or equal to total parts:{} > {} is wrong".format(
                args.part[0], args.part[1]
            )
        )
        sys.exit()
    if args.part[0] < 0:
        logging.critical(
            "Part should be at least zero: {} < 0 is wrong".format(args.part[0])
        )
        sys.exit()

    npzdata = np.load(args.prepfile, encoding="latin1", allow_pickle=True)
    pca_corrected_data = npzdata["pca_corrected_data"]
    masked_bins_per_chr = npzdata["masked_bins_per_chr"]
    masked_bins_per_chr_cum = npzdata["masked_bins_per_chr_cum"]

    indexes, distances, null_ratios = get_reference(
        pca_corrected_data,
        masked_bins_per_chr,
        masked_bins_per_chr_cum,
        ref_size=args.refsize,
        part=args.part[0],
        split_parts=args.part[1],
    )

    np.savez_compressed(
        "{}_{}.npz".format(args.partfile, str(args.part[0])),
        indexes=indexes,
        distances=distances,
        null_ratios=null_ratios,
    )


"""
Merges separate subfiles (one for each thread) to a
new temporary output file.
"""


def tool_newref_post(args, cpus):
    npzdata_prep = np.load(args.prepfile, encoding="latin1", allow_pickle=True)

    big_indexes = []
    big_distances = []
    big_null_ratios = []
    for part in range(1, cpus + 1):
        infile = "{}_{}.npz".format(args.partfile, str(part))
        npzdata_part = np.load(infile, encoding="latin1")
        big_indexes.extend(npzdata_part["indexes"])
        big_distances.extend(npzdata_part["distances"])
        big_null_ratios.extend(npzdata_part["null_ratios"])

    indexes = np.array(big_indexes)
    distances = np.array(big_distances)
    null_ratios = np.array(big_null_ratios)

    np.savez_compressed(
        args.tmpoutfile,
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


def force_remove(file_id):
    attemp = 1
    while True:
        try:
            os.remove(file_id)
            break
        except:
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


def tool_newref_merge(args, outfiles, trained_cutoff):
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
    final_ref["is_nipt"] = args.nipt
    final_ref["trained_cutoff"] = trained_cutoff
    np.savez_compressed(args.outfile, **final_ref)
