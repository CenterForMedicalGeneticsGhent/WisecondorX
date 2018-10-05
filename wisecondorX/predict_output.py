# WisecondorX

import os
import numpy as np


'''
Writes plots.
'''

def exec_write_plots(rem_input, results):
	json_plot_dir = os.path.abspath(rem_input['args'].outid + '_plot_tmp')
	json_dict = {
		'R_script': str('{}/include/plotter.R'.format(rem_input['wd'])),
		'ref_gender': str(rem_input['ref_gender']),
		'beta': str(rem_input['args'].beta),
		'binsize': str(rem_input['binsize']),
		'n_reads': str(rem_input['n_reads']),
		'results_r': results['results_r'],
		'results_w': results['results_w'],
		'results_c': results['results_c'],
		'infile': str('{}.json'.format(json_plot_dir)),
		'out_dir': str('{}.plots'.format(rem_input['args'].outid)),
	}

	from overall_tools import exec_R
	exec_R(json_dict)


'''
Calculates zz-scores, marks aberrations and
writes tables.
'''

def generate_output_tables(rem_input, results):

	_generate_bins_bed(rem_input, results)
	_generate_segments_and_aberrations_bed(rem_input, results)
	_generate_chr_statistics_file(rem_input, results)

def _generate_bins_bed(rem_input, results):
	bins_file = open('{}_bins.bed'.format(rem_input['args'].outid), 'w')
	bins_file.write('chr\tstart\tend\tid\tratio\tzscore\n')
	results_r = results['results_r']
	results_z = results['results_z']
	binsize = rem_input['binsize']

	for chr in range(len(results_r)):
		chr_name = str(chr + 1)
		if chr_name == '23':
			chr_name = 'X'
		if chr_name == '24':
			chr_name = 'Y'
		feat = 1
		for i in range(len(results_r[chr])):
			r = results_r[chr][i]
			z = results_z[chr][i]
			if r == 0:
				r = 'NaN'
			if z == 0:
				z = 'NaN'
			feat_str = '{}:{}-{}'.format(chr_name, str(feat), str(feat + binsize - 1))
			row = [chr_name, feat, feat + binsize - 1, feat_str, r, z]
			bins_file.write('{}\n'.format('\t'.join([str(x) for x in row])))
			feat += binsize
	bins_file.close()

def _generate_segments_and_aberrations_bed(rem_input, results):
	segments_file = open('{}_segments.bed'.format(rem_input['args'].outid), 'w')
	abberations_file = open('{}_aberrations.bed'.format(rem_input['args'].outid), 'w')
	segments_file.write('chr\tstart\tend\tratio\tzscore\n')
	abberations_file.write('chr\tstart\tend\tratio\tzscore\ttype\n')

	for segment in results['results_c']:
		chr_name = str(segment[0] + 1)
		if chr_name == '23':
			chr_name = 'X'
		if chr_name == '24':
			chr_name = 'Y'
		row = [chr_name,
			   int(segment[1] * rem_input['binsize'] + 1),
			   int(segment[2] * rem_input['binsize']),
			   segment[4], segment[3]]
		segments_file.write('{}\n'.format('\t'.join([str(x) for x in row])))
		ploidy = 2

		if (chr_name == 'X' or chr_name == 'Y') and rem_input['actual_gender'] == 'M':
			ploidy = 1
		if float(segment[4]) > __get_aberration_cutoff(rem_input['args'].beta, ploidy)[1]:
			abberations_file.write('{}\tgain\n'.format('\t'.join([str(x) for x in row])))
		elif float(segment[4]) < __get_aberration_cutoff(rem_input['args'].beta, ploidy)[0]:
			abberations_file.write('{}\tloss\n'.format('\t'.join([str(x) for x in row])))

	segments_file.close()
	abberations_file.close()

def __get_aberration_cutoff(beta, ploidy):
	loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
	gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
	return loss_cutoff, gain_cutoff

def _generate_chr_statistics_file(rem_input, results):
	stats_file = open('{}_chr_statistics.txt'.format(rem_input['args'].outid), 'w')
	stats_file.write('chr\tratio.mean\tratio.median\tzscore\n')
	chr_scores = []

	stouffer_scores = []
	for chr in range(len(results['results_z'])):
		stouffer = np.sum(np.array(results['results_z'][chr]) * np.array(results['results_w'][chr])) \
				   / np.sqrt(np.sum(np.power(np.array(results['results_w'][chr]), 2)))
		stouffer_scores.append(stouffer)

	from overall_tools import get_zz_score
	zz_scores = get_zz_score(stouffer_scores)

	for chr in range(len(zz_scores)):

		chr_name = str(chr + 1)
		if chr_name == '23':
			chr_name = 'X'
		if chr_name == '24':
			chr_name = 'Y'
		R = [x for x in results['results_r'][chr] if x != 0]
		chr_ratio_mean = np.mean(R)
		chr_ratio_median = np.median(R)

		row = [chr_name,
			   chr_ratio_median,
			   chr_ratio_mean,
			   zz_scores[chr]]

		stats_file.write('\t'.join([str(x) for x in row]) + '\n')

		chr_scores.append(chr_ratio_mean)

	stats_file.write('Standard deviation ratios (per chromosome): {}\n'
					 .format(str(np.std([x for x in chr_scores if not np.isnan(x)]))))

	stats_file.write('Median segment variance (per bin): {}\n'
					 .format(str(__get_median_within_segment_variance(results['results_c'], results['results_r']))))
	stats_file.close()

def __get_median_within_segment_variance(results_c, results_r):
	vars = []
	for segment in results_c:
		segment_ratios = results_r[segment[0]][int(segment[1]):int(segment[2])]
		segment_ratios = [x for x in segment_ratios if x != 0]
		if segment_ratios:
			var = np.var(segment_ratios)
			vars.append(var)
	return np.median([x for x in vars if not np.isnan(x)])