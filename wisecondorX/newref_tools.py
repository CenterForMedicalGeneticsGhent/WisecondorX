# WisecondorX

import logging
import bisect
import numpy as np
from sklearn.decomposition import PCA


'''
Finds mask (locations of bins without data) in the
subset 'samples'.
'''

def get_mask(samples):
	by_chr= []
	bins_per_chr = []
	sample_count = len(samples)

	for chr in range(1, 25):
		max_len = max([sample[str(chr)].shape[0] for sample in samples])
		this_chr= np.zeros((max_len, sample_count), dtype=float)
		bins_per_chr.append(max_len)
		i = 0
		for sample in samples:
			this_chr[:, i] = sample[str(chr)]
			i += 1
		by_chr.append(this_chr)
	all_data = np.concatenate(by_chr, axis=0)

	sum_per_sample = np.sum(all_data, 0)
	all_data = all_data / sum_per_sample

	sum_per_bin = np.sum(all_data, 1)
	mask = sum_per_bin > 0

	return mask, bins_per_chr


'''
Normalizes samples for read depth and applies mask.
'''

def normalize_and_mask(samples, chrs, mask):
	by_chr= []
	sample_count = len(samples)

	for chr in chrs:
		max_len = max([sample[str(chr)].shape[0] for sample in samples])
		this_chr= np.zeros((max_len, sample_count), dtype=float)
		i = 0
		for sample in samples:
			this_chr[:, i] = sample[str(chr)]
			i += 1
		by_chr.append(this_chr)
	all_data = np.concatenate(by_chr, axis=0)

	sum_per_sample = np.sum(all_data, 0)
	all_data = all_data / sum_per_sample

	masked_data = all_data[mask, :]

	return masked_data


'''
Executes PCA. Rotations are saved which enable
between sample normalization in the test phase.
'''

def train_pca(ref_data, pcacomp=3):
	t_data = ref_data.T
	pca = PCA(n_components=pcacomp)
	pca.fit(t_data)
	PCA(copy=True, whiten=False)
	transformed = pca.transform(t_data)
	inversed = pca.inverse_transform(transformed)
	corrected = t_data / inversed

	return corrected.T, pca


'''
Calculates within-sample reference.
'''

def get_reference(pca_corrected_data, masked_bins_per_chr, masked_bins_per_chr_cum,
				   ref_size, part, split_parts):
	big_indexes = []
	big_distances = []

	bincount = masked_bins_per_chr_cum[-1]

	start_num, end_num = _get_part(part - 1, split_parts, bincount)
	logging.info('Working on thread {} of {}, meaning bins {} up to {}'
				 .format(part,split_parts, start_num, end_num))
	regions = _split_by_chr(start_num, end_num, masked_bins_per_chr_cum)

	for region in regions:
		chr= region[0]
		start = region[1]
		end = region[2]

		if start_num > start:
			start = start_num
		if end_num < end:
			end = end_num

		if len(masked_bins_per_chr_cum) > 22 and chr!= 22 and chr!= 23:
			part_indexes = np.zeros((end - start, ref_size), dtype=np.int32)
			part_distances = np.ones((end - start, ref_size))
			big_indexes.extend(part_indexes)
			big_distances.extend(part_distances)
			continue
		chr_data = np.concatenate((pca_corrected_data[:masked_bins_per_chr_cum[chr] -
													   masked_bins_per_chr[chr], :],
									 pca_corrected_data[masked_bins_per_chr_cum[chr]:, :]))

		part_indexes, part_distances = get_ref_for_bins(ref_size, start, end,
														 pca_corrected_data, chr_data)

		big_indexes.extend(part_indexes)
		big_distances.extend(part_distances)

	index_array = np.array(big_indexes)
	distance_array = np.array(big_distances)

	return index_array, distance_array

def _split_by_chr(start, end, chr_bin_sums):
	areas = []
	tmp = [0, start, 0]
	for i, val in enumerate(chr_bin_sums):
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

def _get_part(partnum, outof, bincount):
	start_bin = int(bincount / float(outof) * partnum)
	end_bin = int(bincount / float(outof) * (partnum + 1))
	return start_bin, end_bin


'''
Calculates within-sample reference for a particular chromosome.
'''

def get_ref_for_bins(ref_size, start, end, pca_corrected_data, chr_data):
	find_pos = bisect.bisect
	ref_indexes = np.zeros((end - start, ref_size), dtype=np.int32)
	ref_distances = np.ones((end - start, ref_size))
	for this_bin in range(start, end):
		this_mask = np.sum(np.power(chr_data - pca_corrected_data[this_bin, :], 2), 1)
		this_indexes = [-1 for i in range(ref_size)]
		this_distances = [1e10 for i in range(ref_size)]
		remove_index = this_indexes.pop
		remove_dist = this_distances.pop
		insert_index = this_indexes.insert
		insert_dist = this_distances.insert
		cur_max = 1e10
		for i, binVal in enumerate(this_mask):
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