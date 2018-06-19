import bisect
import json
import logging
import math
import os
import sys
import subprocess
import time
import pysam
import numpy as np
from sklearn.decomposition import PCA
import warnings

warnings.filterwarnings('ignore', 'Mean of empty slice')
warnings.filterwarnings('ignore', 'Degrees of freedom <= 0 for slice')

np.seterr('ignore')
np_sum = np.sum
np_pow = np.power
np_max = np.argmax
np_mean = np.mean
np_median = np.median
np_std = np.std
np_abs = np.abs
np_sqrt = np.sqrt

find_pos = bisect.bisect


def train_pca(ref_data, pcacomp=3):
    t_data = ref_data.T
    pca = PCA(n_components=pcacomp)
    pca.fit(t_data)
    PCA(copy=True, whiten=False)
    transformed = pca.transform(t_data)
    inversed = pca.inverse_transform(transformed)
    corrected = t_data / inversed

    return corrected.T, pca


def apply_pca(sample_data, mean, components):
    pca = PCA(n_components=components.shape[0])
    pca.components_ = components
    pca.mean_ = mean

    transform = pca.transform(np.array([sample_data]))

    reconstructed = np.dot(transform, pca.components_) + pca.mean_
    reconstructed = reconstructed[0]
    return sample_data / reconstructed


def convert_bam(bamfile, binsize, min_shift, threshold, mapq=1, demand_pair=False):
    # Prepare the list of chromosomes
    chromosomes = dict()
    for chromosome in range(1, 25):
        chromosomes[str(chromosome)] = None

    # Flush the current stack of reads
    def flush(read_buff, counts):
        stair_size = len(read_buff)
        if stair_size <= threshold or threshold < 0:
            for read in read_buff:
                location = read.pos / binsize
                counts[int(location)] += 1

    sam_file = pysam.AlignmentFile(bamfile, "rb")
    reads_seen = 0
    reads_kept = 0
    reads_mapq = 0
    reads_rmdup = 0
    reads_pairf = 0
    larp = -1  # Last Read Position...
    larp2 = -1

    for index, chrom in enumerate(sam_file.references):

        chrom_name = chrom
        if chrom_name[:3].lower() == 'chr':
            chrom_name = chrom_name[3:]
        if chrom_name not in chromosomes and chrom_name != "X" and chrom_name != "Y":
            continue

        logging.info('Working at {}; processing {} bins'
                     .format(chrom, int(sam_file.lengths[index] / float(binsize) + 1)))
        counts = np.zeros(int(sam_file.lengths[index] / float(binsize) + 1), dtype=np.int32)

        read_buff = []
        sam_iter = sam_file.fetch(chrom)

        if chrom_name == 'X':
            chrom_name = '23'
        if chrom_name == 'Y':
            chrom_name = '24'
        prev_read = sam_iter.next()
        # Split paths here, for-loop was heavily slowed down by if-statements otherwise
        if demand_pair:
            for read in sam_iter:
                if (int(read.pos) - int(prev_read.pos)) > min_shift:
                    flush(read_buff, counts)
                    read_buff = []
                # Normal ndups will be appended here

                if not (read.is_proper_pair and read.is_read1):
                    reads_pairf += 1
                    continue

                if larp == read.pos and larp2 == read.next_reference_start:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= mapq:
                        read_buff.append(read)
                        prev_read = read
                    else:
                        reads_mapq += 1

                larp2 = read.next_reference_start

                reads_seen += 1
                larp = read.pos
        else:
            for read in sam_iter:
                if (int(read.pos) - int(prev_read.pos)) > min_shift:
                    flush(read_buff, counts)
                    read_buff = []
                # Normal ndups will be appended here

                if larp == read.pos:
                    reads_rmdup += 1
                else:
                    if read.mapping_quality >= mapq:
                        read_buff.append(read)
                        prev_read = read
                    else:
                        reads_mapq += 1

                reads_seen += 1
                larp = read.pos

            # Flush after we're done
            flush(read_buff, counts)
            chromosomes[chrom_name] = counts
            reads_kept += sum(counts)

    # print reads_seen,reads_kept
    qual_info = {'mapped': sam_file.mapped,
                 'unmapped': sam_file.unmapped,
                 'no_coordinate': sam_file.nocoordinate,
                 'filter_rmdup': reads_rmdup,
                 'filter_mapq': reads_mapq,
                 'pre_retro': reads_seen,
                 'post_retro': reads_kept,
                 'pair_fail': reads_pairf}
    return chromosomes, qual_info


def scale_sample(sample, from_size, to_size):
    if not to_size or from_size == to_size:
        return sample

    if to_size == 0 or from_size == 0 or to_size < from_size or to_size % from_size > 0:
        logging.error("Impossible binsize scaling requested: {} to {}".format(int(from_size), int(to_size)))
        sys.exit()

    return_sample = dict()
    scale = to_size / from_size
    for chrom in sample:
        chrom_data = sample[chrom]
        new_len = int(np.ceil(len(chrom_data) / float(scale)))
        scaled_chrom = np.zeros(new_len, dtype=np.int32)
        for i in range(new_len):
            scaled_chrom[i] = np_sum(chrom_data[int(i * scale):int(i * scale + scale)])
            return_sample[chrom] = scaled_chrom
    return return_sample


def to_numpy_array(samples, gender):
    by_chrom = []
    chrom_bins = []
    sample_count = len(samples)

    if gender == "F":
        chromosomes = range(1, 24)
    else:
        chromosomes = range(1, 25)

    for chromosome in chromosomes:
        max_len = max([sample[str(chromosome)].shape[0] for sample in samples])
        this_chrom = np.zeros((max_len, sample_count), dtype=float)
        chrom_bins.append(max_len)
        i = 0
        for sample in samples:
            this_chrom[:, i] = sample[str(chromosome)]
            i += 1
        by_chrom.append(this_chrom)
    all_data = np.concatenate(by_chrom, axis=0)

    sum_per_sample = np_sum(all_data, 0)
    all_data = all_data / sum_per_sample

    logging.info('Applying nonzero mask on the data: {}'.format(all_data.shape))
    sum_per_bin = np_sum(all_data, 1)
    mask = sum_per_bin > 0
    masked_data = all_data[mask, :]
    logging.info('becomes {}'.format(masked_data.shape))

    return masked_data, chrom_bins, mask


def to_numpy_ref_format(sample, chrom_bins, mask, gender):
    by_chrom = []

    if gender == "M":
        chrs = range(1, 25)
    else:
        chrs = range(1, 24)

    for chromosome in chrs:
        this_chrom = np.zeros(chrom_bins[chromosome - 1], dtype=float)
        min_len = min(chrom_bins[chromosome - 1], len(sample[str(chromosome)]))
        this_chrom[:min_len] = sample[str(chromosome)][:min_len]
        by_chrom.append(this_chrom)
    all_data = np.concatenate(by_chrom, axis=0)
    all_data = all_data / np_sum(all_data)
    masked_data = all_data[mask]

    return masked_data


def inflate_array(array, mask):
    temp = np.zeros(mask.shape[0])
    j = 0
    for i, val in enumerate(mask):
        if val:
            temp[i] = array[j]
            j += 1
    return temp


def inflate_array_multi(array, mask_list):
    temp = array
    for mask in reversed(mask_list):
        temp = inflate_array(temp, mask)
    return temp


def get_ref_for_bins(amount, start, end, sample_data, other_data, knit_length):
    ref_indexes = np.zeros((end - start, amount), dtype=np.int32)
    ref_distances = np.ones((end - start, amount))
    for this_bin in xrange(start, end):
        this_mask = np_sum(np_pow(other_data - sample_data[this_bin, :], 2), 1)
        this_indexes = [-1 for i in xrange(amount)]
        this_distances = [1e10 for i in xrange(amount)]
        remove_index = this_indexes.pop
        remove_dist = this_distances.pop
        insert_index = this_indexes.insert
        insert_dist = this_distances.insert
        cur_max = 1e10
        tot_len = len(this_mask) + (end - start)
        for i, binVal in enumerate(this_mask):
            if (end - start) + i >= tot_len - knit_length:  # skip XY referenced bins
                continue
            if binVal < cur_max:
                pos = find_pos(this_distances, binVal)
                remove_index(-1)
                remove_dist(-1)
                insert_index(pos, i)
                insert_dist(pos, binVal)
                cur_max = this_distances[-1]
        ref_indexes[this_bin - start, :] = this_indexes
        ref_distances[this_bin - start, :] = this_distances
    return ref_indexes, ref_distances


def get_optimal_cutoff(reference, chromosome_sizes, repeats):
    autosome_len = sum(chromosome_sizes[:22])

    autosome_ref = reference[:autosome_len]
    autosome_cutoff = float("inf")
    for i in range(0, repeats):
        mask = autosome_ref < autosome_cutoff
        average = np.average(autosome_ref[mask])
        stddev = np.std(autosome_ref[mask])
        autosome_cutoff = average + 3 * stddev

    allosome_ref = reference[autosome_len:]
    allosome_cutoff = float("inf")
    for i in range(0, repeats):
        mask = allosome_ref < allosome_cutoff
        average = np.average(allosome_ref[mask])
        stddev = np.std(allosome_ref[mask])
        allosome_cutoff = average + 3 * stddev

    return autosome_cutoff, allosome_cutoff


# Returns: Chromosome index, startBinNumber, endBinNumber
def split_by_chrom(start, end, chromosome_bin_sums):
    areas = []
    tmp = [0, start, 0]
    for i, val in enumerate(chromosome_bin_sums):
        tmp[0] = i
        if val >= end:
            break
        if start < val < end:
            tmp[2] = val
            areas.append(tmp)
            tmp = [i, val, 0]
        tmp[1] = val
    tmp[2] = end
    areas.append(tmp)
    return areas


# Returns: Start and end bin numbers this instance should work on
def get_part(partnum, outof, bincount):
    start_bin = int(bincount / float(outof) * partnum)
    end_bin = int(bincount / float(outof) * (partnum + 1))
    return start_bin, end_bin


def get_reference(corrected_data, chromosome_bins, chromosome_bin_sums,
                  gender, select_ref_amount=100, part=1, split_parts=1):
    time_start_selection = time.time()
    big_indexes = []
    big_distances = []

    bincount = chromosome_bin_sums[-1]

    start_num, end_num = get_part(part - 1, split_parts, bincount)
    logging.info('Working on thread {} of {} meaning bins {} up to {}'.format(part, split_parts, start_num, end_num))
    regions = split_by_chrom(start_num, end_num, chromosome_bin_sums)

    for region in regions:
        chrom = region[0]
        start = region[1]
        end = region[2]

        if start_num > start:
            start = start_num
        if end_num < end:
            end = end_num

        logging.info('Thread {} | Actual Chromosome area {} {} | chr {}'.format(
            part, chromosome_bin_sums[chrom] - chromosome_bins[chrom], chromosome_bin_sums[chrom], str(chrom + 1)))
        chrom_data = np.concatenate((corrected_data[:chromosome_bin_sums[chrom] - chromosome_bins[chrom], :],
                                    corrected_data[chromosome_bin_sums[chrom]:, :]))

        x_length = chromosome_bin_sums[22] - (chromosome_bin_sums[22] - chromosome_bins[22])  # index 22 -> chrX

        if gender == "M":
            y_length = chromosome_bin_sums[23] - (chromosome_bin_sums[23] - chromosome_bins[23])
            if chrom == 22:
                knit_length = y_length
            elif chrom == 23:
                knit_length = x_length
            else:
                knit_length = x_length + y_length
        else:
            if chrom == 22:
                knit_length = 0
            else:
                knit_length = x_length

        part_indexes, part_distances = get_ref_for_bins(select_ref_amount, start,
                                                        end, corrected_data, chrom_data, knit_length)
        big_indexes.extend(part_indexes)
        big_distances.extend(part_distances)

        logging.info('{} Time spent {} seconds'.format(part, int(time.time() - time_start_selection)))

    index_array = np.array(big_indexes)
    distance_array = np.array(big_distances)

    return index_array, distance_array


def try_sample(test_data, test_copy, indexes, distances, chromosome_bins,
               chromosome_bin_sums, autosome_cutoff, allosome_cutoff):
    bincount = chromosome_bin_sums[-1]

    results_z = np.zeros(bincount)
    results_r = np.zeros(bincount)
    ref_sizes = np.zeros(bincount)
    std_dev_sum = 0.
    std_dev_num = 0
    i = 0

    for chrom in xrange(len(chromosome_bins)):
        start = chromosome_bin_sums[chrom] - chromosome_bins[chrom]
        end = chromosome_bin_sums[chrom]
        chrom_data = np.concatenate(
            (test_copy[:chromosome_bin_sums[chrom] - chromosome_bins[chrom]], test_copy[chromosome_bin_sums[chrom]:]))

        for index in indexes[start:end]:
            if chrom < 22:
                ref_data = chrom_data[index[distances[i] < autosome_cutoff]]
            else:
                ref_data = chrom_data[index[distances[i] < allosome_cutoff]]
            ref_data = ref_data[ref_data >= 0]  # Previously found aberrations are marked by negative values
            ref_mean = np_mean(ref_data)
            ref_stdev = np_std(ref_data)
            if not np.isnan(ref_stdev):
                std_dev_sum += ref_stdev
                std_dev_num += 1
            results_z[i] = (test_data[i] - ref_mean) / ref_stdev
            results_r[i] = test_data[i] / ref_mean
            ref_sizes[i] = ref_data.shape[0]
            i += 1

    return results_z, results_r, ref_sizes, std_dev_sum / std_dev_num


def repeat_test(test_data, indexes, distances, chromosome_bins,
                chromosome_bin_sums, autosome_cutoff, allosome_cutoff, threshold, repeats):
    results_z = None
    results_r = None
    test_copy = np.copy(test_data)
    for i in xrange(repeats):
        results_z, results_r, ref_sizes, std_dev_avg = try_sample(test_data, test_copy, indexes, distances,
                                                                  chromosome_bins, chromosome_bin_sums,
                                                                  autosome_cutoff, allosome_cutoff)
        test_copy[np_abs(results_z) >= threshold] = -1
    return results_z, results_r, ref_sizes, std_dev_avg


def get_weights(distances):
    inverse_weights = [np.mean(x) for x in distances]
    weights = np.array([1/x for x in inverse_weights])
    return weights


def generate_txt_output(args, json_out):
    bed_file = open(args.outid + "_bins.bed", "w")
    bed_file.write("chr\tstart\tend\tid\tzscore\tratio\n")
    results_r = json_out["results_r"]
    results_z = json_out["results_z"]
    results_w = json_out["results_w"]
    binsize = json_out["binsize"]
    for chr_i in range(len(results_r)):
        chr = str(chr_i + 1)
        if chr == "23":
            chr = "X"
        if chr == "24":
            chr = "Y"
        feat = 1
        for feat_i in range(len(results_r[chr_i])):
            r = results_r[chr_i][feat_i]
            z = results_z[chr_i][feat_i]
            if r == 0:
                r = "NaN"
            if z == 0:
                z = "NaN"
            feat_str = chr + ":" + str(feat) + "-" + str(feat + binsize - 1)
            it = [chr, feat, feat + binsize - 1, feat_str, z, r]
            it = [str(x) for x in it]
            bed_file.write("\t".join(it) + "\n")
            feat += binsize
    bed_file.close()

    segments_file = open(args.outid + "_segments.bed", "w")
    ab_file = open(args.outid + "_aberrations.bed", "w")
    segments_file.write("chr\tstart\tend\tzscore\tratio\n")
    ab_file.write("chr\tstart\tend\tzscore\tratio\ttype\n")
    segments = json_out["cbs_calls"]
    for segment in segments:
        chr = str(int(segment[0]))
        if chr == "23":
            chr = "X"
        if chr == "24":
            chr = "Y"
        it = [chr, int(segment[1] * binsize + 1), int((segment[2] + 1) * binsize), segment[3], segment[4]]
        it = [str(x) for x in it]
        segments_file.write("\t".join(it) + "\n")
        if float(segment[4]) > np.log2(1. + args.beta / 4):
            ab_file.write("\t".join(it) + "\tgain\n")
        elif float(segment[4]) < np.log2(1. - args.beta / 4):
            ab_file.write("\t".join(it) + "\tloss\n")

    segments_file.close()

    statistics_file = open(args.outid + "_chr_statistics.txt", "w")
    statistics_file.write("chr\tzscore\tratio.mean\tratio.median\n")
    chrom_scores = []
    for chr_i in range(len(results_r)):
        chr = str(chr_i + 1)
        if chr == "23":
            chr = "X"
        if chr == "24":
            chr = "Y"
        R = [x for x in results_r[chr_i] if x != 0]

        stouffer = np.sum(np.array(results_z[chr_i]) * np.array(results_w[chr_i])) \
                   / np.sqrt(np.sum(np.power(np.array(results_w[chr_i]), 2)))

        chrom_ratio_mean = np.mean(R)
        chrom_ratio_median = np.median(R)

        statistics_file.write(str(chr)
                              + "\t" + str(stouffer)
                              + "\t" + str(chrom_ratio_median)
                              + "\t" + str(chrom_ratio_mean)
                              + "\n")
        chrom_scores.append(chrom_ratio_mean)

    statistics_file.write("Standard deviation mean chromosomal ratio: " +
                          str(np.std([x for x in chrom_scores if not np.isnan(x)])) + "\n")
    statistics_file.write("Median within-segment binwise ratio variance: " + str(
        get_median_within_segment_variance(segments, results_r)) + "\n")
    statistics_file.close()


def write_plots(args, json_out, wc_dir, gender):
    json_file = open(args.outid + "_plot_tmp.json", "w")
    json.dump(json_out, json_file)
    json_file.close()

    plot_script = str(os.path.dirname(wc_dir)) + "/include/plotter.R"
    if gender == "M":
        sexchrom = "XY"
    else:
        sexchrom = "X"
    r_cmd = ["Rscript", plot_script,
             "--infile", "{}_plot_tmp.json".format(args.outid),
             "--outdir", args.outid + ".plots",
             "--sexchroms", sexchrom,
             "--beta", str(args.beta)]
    logging.debug("plot cmd: {}".format(r_cmd))

    try:
        subprocess.check_call(r_cmd)
    except subprocess.CalledProcessError as e:
        logging.critical("Script {} failed with error {}".format(plot_script, e))
        sys.exit()

    os.remove(args.outid + "_plot_tmp.json")


def get_median_within_segment_variance(segments, binratios):
    vars = []
    for segment in segments:
        segment_ratios = binratios[int(segment[0]) - 1][int(segment[1]):int(segment[2])]
        segment_ratios = [x for x in segment_ratios if x != 0]
        if [] != segment_ratios:
            var = np.var(segment_ratios)
            vars.append(var)
    return np.median([x for x in vars if not np.isnan(x)])


def apply_blacklist(args, binsize, results_r, results_z, results_w, sample, gender):
    blacklist = {}

    for line in open(args.blacklist):
        bchr, bstart, bstop = line.strip().split("\t")
        bchr = bchr[3:]
        if bchr not in blacklist.keys():
            blacklist[bchr] = []
        blacklist[bchr].append([int(int(bstart) / binsize), int(int(bstop) / binsize) + 1])

    for chr in blacklist.keys():
        for s_s in blacklist[chr]:
            if chr == "X":
                chr = "23"
            if chr == "Y":
                chr = "24"
            for pos in range(s_s[0], s_s[1]):
                if gender == "F" and chr == "24":
                    continue
                results_r[int(chr) - 1][pos] = 0
                results_z[int(chr) - 1][pos] = 0
                results_w[int(chr) - 1][pos] = 0
                sample[chr][pos] = 0


def cbs(args, results_r, results_z, results_w, gender, wc_dir):
    json_cbs_temp_dir = os.path.abspath(args.outid + "_CBS_tmp")
    json_cbs_file = open(json_cbs_temp_dir + "_01.json", "w")
    json.dump({"results_r": results_r, "weights": results_w}, json_cbs_file)
    json_cbs_file.close()
    cbs_script = str(os.path.dirname(wc_dir)) + "/include/CBS.R"

    if gender == "M":
        sexchrom = "XY"
    else:
        sexchrom = "X"
    r_cmd = ["Rscript", cbs_script,
             "--infile", "{}_01.json".format(json_cbs_temp_dir),
             "--outfile", "{}_02.json".format(json_cbs_temp_dir),
             "--sexchroms", sexchrom,
             "--alpha", str(args.alpha)]
    logging.debug("CBS cmd: {}".format(r_cmd))

    try:
        subprocess.check_call(r_cmd)
    except subprocess.CalledProcessError as e:
        logging.critical("Script {} failed with error {}".format(cbs_script, e))
        sys.exit()
    os.remove(json_cbs_temp_dir + "_01.json")
    cbs_data = json.load(open(json_cbs_temp_dir + "_02.json"))[1:]
    cbs_data = [[float(y.encode("utf-8")) for y in x] for x in cbs_data]
    os.remove(json_cbs_temp_dir + "_02.json")

    stouffer_scores = []
    for cbs_call_index in range(len(cbs_data[0])):
        chr_i = int(cbs_data[0][cbs_call_index]) - 1
        start = int(cbs_data[1][cbs_call_index]) - 1
        end = int(cbs_data[2][cbs_call_index])  # no - 1! (closed interval in python)

        z_segment = np.array(results_z[chr_i][start:end])
        w_segment = np.array(results_w[chr_i][start:end])

        stouffer = np.sum(z_segment * w_segment) / np.sqrt(np.sum(np.power(w_segment,2)))
        stouffer_scores.append(stouffer)

    # Save results

    cbs_calls = []
    for cbs_call_index in range(len(cbs_data[0])):
        cbs_calls.append(
            [cbs_data[0][cbs_call_index], cbs_data[1][cbs_call_index] - 1, cbs_data[2][cbs_call_index] - 1,
             stouffer_scores[cbs_call_index], cbs_data[4][cbs_call_index]])
    return cbs_calls
