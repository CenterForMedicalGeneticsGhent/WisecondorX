#!/usr/bin/env python3

import argparse
import logging
import warnings
import sys
import os
import numpy as np


def tool_convert(args):

	logging.info('Starting conversion')

	from convert_tools import convert_bam
	sample, qual_info = convert_bam(args)

	if args.gender:
		gender = args.gender
		logging.info('Gender {}'.format(gender))
	else:
		logging.info('Predicting gender ...')
		from convert_tools import get_gender
		gender = get_gender(args, sample)
		logging.info('... {}'.format(gender))

	np.savez_compressed(args.outfile,
						binsize=args.binsize,
						sample=sample,
						gender=gender,
						quality=qual_info)

	logging.info('Finished conversion')


def tool_newref(args):

	logging.info('Creating new reference')

	split_path = list(os.path.split(args.outfile))
	if split_path[-1][-4:] == '.npz':
		split_path[-1] = split_path[-1][:-4]
	base_path = os.path.join(split_path[0], split_path[1])

	args.basepath = base_path
	args.prepfile = '{}_prep.npz'.format(base_path)
	args.partfile = '{}_part'.format(base_path)

	samples = []
	genders = []
	logging.info('Importing data ...')
	for infile in args.infiles:
		logging.info('Loading: {}'.format(infile))
		npzdata = np.load(infile, encoding='latin1')
		sample = npzdata['sample'].item()
		binsize = int(npzdata['binsize'])
		gender = str(npzdata['gender'], 'utf-8').rstrip('\x00')
		logging.info('Binsize: {} | Gender : {}'.format(int(binsize), gender))

		from overall_tools import scale_sample, gender_correct
		scaled_sample = scale_sample(sample, binsize, args.binsize)
		corrected_sample = gender_correct(scaled_sample, gender)

		genders.append(gender)
		samples.append(corrected_sample)

	samples = np.array(samples)

	from newref_tools import get_mask
	total_mask, bins_per_chr = get_mask(samples)
	if genders.count('F') > 4:
		mask_F, _ = get_mask(samples[np.array(genders) == 'F'])
		total_mask = total_mask & mask_F
	if genders.count('M') > 4:
		mask_M, _ = get_mask(samples[np.array(genders) == 'M'])
		total_mask = total_mask & mask_M

	from newref_control import tool_newref_prep, tool_newref_main, tool_newref_merge
	outfiles = []
	if len(genders) > 9:
		logging.info('Starting autosomal reference creation ...')
		args.tmpoutfile = '{}.tmp.A.npz'.format(args.basepath)
		outfiles.append(args.tmpoutfile)
		tool_newref_prep(args, samples, 'A', total_mask, bins_per_chr)
		logging.info('This might take a while ...')
		tool_newref_main(args, args.cpus)
	else:
		logging.critical('Provide at least 10 samples to enable the generation of a reference.')
		sys.exit()

	if genders.count('F') > 4:
		logging.info('Starting female gonosomal reference creation ...')
		args.tmpoutfile = '{}.tmp.F.npz'.format(args.basepath)
		outfiles.append(args.tmpoutfile)
		tool_newref_prep(args, samples[np.array(genders) == 'F'], 'F', total_mask, bins_per_chr)
		logging.info('This might take a while ...')
		tool_newref_main(args, 1)
	else:
		logging.warning('Provide at least 5 female samples to enable normalization of female gonosomes.')

	if genders.count('M') > 4:
		logging.info('Starting male gonosomal reference creation ...')
		args.tmpoutfile = '{}.tmp.M.npz'.format(args.basepath)
		outfiles.append(args.tmpoutfile)
		tool_newref_prep(args, samples[np.array(genders) == 'M'], 'M', total_mask, bins_per_chr)
		tool_newref_main(args, 1)
	else:
		logging.warning('Provide at least 5 male samples to enable normalization of male gonosomes. '
						'If these are of no interest (e.g. NIPT), ignore this warning.')

	tool_newref_merge(args, outfiles)

	logging.info('Finished creating reference')


def tool_test(args):

	logging.info('Starting CNA prediction')

	if not args.bed and not args.plot:
		logging.critical('No output format selected. '
						 'Select at least one of the supported output formats (--bed, --plot)')
		sys.exit()

	if args.beta <= 0 or args.beta > 1:
		logging.critical('Parameter beta should be a strictly positive number lower than 1')
		sys.exit()

	if args.alpha <= 0 or args.alpha > 1:
		logging.critical('Parameter alpha should be a strictly positive number lower than 1')
		sys.exit()

	if args.beta < 0.05:
		logging.warning('Parameter beta seems to be a bit low. '
						'Have you read https://github.com/CenterForMedicalGeneticsGhent/wisecondorX#parameters '
						'on parameter optimization?')

	if args.alpha > 0.1:
		logging.warning('Parameter alpha seems to be a bit high. '
						'Have you read https://github.com/CenterForMedicalGeneticsGhent/wisecondorX#parameters '
						'on parameter optimization?')

	logging.info('Importing data ...')
	ref_file = np.load(args.reference, encoding='latin1')
	sample_file = np.load(args.infile, encoding='latin1')

	sample = sample_file['sample'].item()
	n_reads = sum([sum(sample[x]) for x in sample.keys()])
	actual_gender = str(sample_file['gender'], 'utf-8').rstrip('\x00')
	ref_gender = actual_gender

	from overall_tools import scale_sample, gender_correct
	sample_bin_size = int(sample_file['binsize'].item())
	sample = scale_sample(sample, sample_bin_size, int(ref_file['binsize']))
	sample = gender_correct(sample, actual_gender)

	logging.info('Normalizing autosomes ...')

	from predict_control import normalize
	results_r, results_z, results_w, ref_sizes = normalize(args, sample, ref_file, 'A')

	if not ref_file['has_male'] and actual_gender == 'M':
		logging.warning('This sample is male, whilst the reference is created with fewer than 5 males. '
						'The female gonosomal reference will be used for X predictions. Note that these might '
						'not be accurate. If the latter is desired, create a new reference and include more '
						'male samples.')
		ref_gender = 'F'

	elif not ref_file['has_female'] and actual_gender == 'F':
		logging.warning('This sample is female, whilst the reference is created with fewer than 5 females. '
						'The male gonosomal reference will be used for XY predictions. Note that these might '
						'not be accurate. If the latter is desired, create a new reference and include more '
						'female samples.')
		ref_gender = 'M'

	logging.info('Normalizing gonosomes ...')

	results_r_2, results_z_2, results_w_2, ref_sizes_2 = normalize(args, sample, ref_file, ref_gender)

	rem_input = {'args': args,
				 'wd' : str(os.path.dirname(os.path.realpath(__file__))),
				 'binsize' : int(ref_file['binsize']),
				 'n_reads': n_reads,
				 'ref_gender': ref_gender,
				 'actual_gender': actual_gender,
				 'mask' : ref_file['mask.{}'.format(ref_gender)],
				 'bins_per_chr' : ref_file['bins_per_chr.{}'.format(ref_gender)],
				 'masked_bins_per_chr' : ref_file['masked_bins_per_chr.{}'.format(ref_gender)],
				 'masked_bins_per_chr_cum' : ref_file['masked_bins_per_chr_cum.{}'.format(ref_gender)]}

	del ref_file

	results_r = np.append(results_r, results_r_2)
	results_z = np.append(results_z, results_z_2)
	results_w = np.append(results_w * np.nanmedian(results_w_2), results_w_2 * np.nanmedian(results_w))
	results_w = results_w / np.nanmedian(results_w)
	ref_sizes = np.append(ref_sizes, ref_sizes_2)

	results = {'results_r': results_r,
			   'results_z': results_z,
			   'results_w': results_w}

	from predict_control import get_post_processed_result
	for result in results.keys():
		results[result] = get_post_processed_result(args, results[result], ref_sizes, rem_input)

	from predict_tools import log_trans
	log_trans(results)

	if args.blacklist:
		logging.info('Applying blacklist ...')
		from predict_tools import apply_blacklist
		apply_blacklist(rem_input, results)

	logging.info('Executing circular binary segmentation ...')

	from predict_tools import exec_cbs
	results['results_c'] = exec_cbs(rem_input, results)

	if args.bed:
		logging.info('Writing tables ...')
		from predict_output import generate_output_tables
		generate_output_tables(rem_input, results)

	if args.plot:
		logging.info('Writing plots ...')
		from predict_output import exec_write_plots
		exec_write_plots(rem_input, results)

	logging.info('Finished prediction')


def output_gender(args):

	npzdata = np.load(args.infile, encoding='latin1')
	gender = str(npzdata['gender'])
	if gender == 'M':
		print('male')
	else:
		print('female')


def main():
	parser = argparse.ArgumentParser(description='wisecondorX')
	parser.add_argument('--loglevel',
						type=str,
						default='INFO',
						choices=['info', 'warning', 'debug', 'error', 'critical'])
	subparsers = parser.add_subparsers()

	parser_convert = subparsers.add_parser('convert',
										   description='Convert and filter a .bam file to a .npz',
										   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_convert.add_argument('infile',
								type=str,
								help='.bam input file for conversion')
	parser_convert.add_argument('outfile',
								type=str,
								help='Output .npz file')
	parser_convert.add_argument('--binsize',
								type=float,
								default=5e3,
								help='Bin size (bp)')
	parser_convert.add_argument('--retdist',
								type=int,
								default=4,
								help='Maximum amount of base pairs difference between sequential reads '
									 'to consider them part of the same tower')
	parser_convert.add_argument('--retthres',
								type=int,
								default=4,
								help='Threshold for when a group of reads is considered a tower and will be removed')
	parser_convert.add_argument('--gender',
								type=str,
								choices=['F', 'M'],
								help='Gender of the case. If not given, WisecondorX will predict it')
	parser_convert.add_argument('--paired',
								action='store_true',
								help='Use paired-end reads | default is single-end')
	parser_convert.add_argument('--ycutoff',
								type=float,
								default=0.004,
								help='A cutoff value representing the ratio \'Y reads/total reads\'. '
									 'Used to predict gender. Might require training for selecting the optimal value')
	parser_convert.set_defaults(func=tool_convert)

	parser_gender = subparsers.add_parser('gender',
										  description='Returns the gender of a .npz resulting from convert',
										  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_gender.add_argument('infile',
							   type=str,
							   help='.npz input file')
	parser_gender.set_defaults(func=output_gender)

	parser_newref = subparsers.add_parser('newref',
										  description='Create a new reference using healthy reference samples',
										  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_newref.add_argument('infiles',
							   type=str,
							   nargs='+',
							   help='Path to all reference data files (e.g. path/to/reference/*.npz)')
	parser_newref.add_argument('outfile',
							   type=str,
							   help='Path and filename for the reference output (e.g. path/to/myref.npz)')
	parser_newref.add_argument('--refsize',
							   type=int,
							   default=300,
							   help='Amount of reference locations per target')
	parser_newref.add_argument('--binsize',
							   type=int,
							   default=1e5,
							   help='Scale samples to this binsize, multiples of existing binsize only')
	parser_newref.add_argument('--cpus',
							   type=int,
							   default=1,
							   help='Use multiple cores to find reference bins')
	parser_newref.set_defaults(func=tool_newref)

	parser_test = subparsers.add_parser('predict',
										description='Find copy number aberrations',
										formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_test.add_argument('infile',
							 type=str,
							 help='.npz input file')
	parser_test.add_argument('reference',
							 type=str,
							 help='Reference .npz, as previously created with newref')
	parser_test.add_argument('outid',
							 type=str,
							 help='Basename (w/o extension) of output files (paths are allowed, e.g. path/to/ID_1)')
	parser_test.add_argument('--minrefbins',
							 type=int,
							 default=150,
							 help='Minimum amount of sensible reference bins per target bin.')
	parser_test.add_argument('--maskrepeats',
							 type=int,
							 default=5,
							 help='Regions with distances > mean + sd * 3 will be masked. Number of masking cycles.')
	parser_test.add_argument('--alpha',
							 type=float,
							 default=1e-4,
							 help='p-value cut-off for calling a CBS breakpoint.')
	parser_test.add_argument('--beta',
							 type=float,
							 default=0.1,
							 help='Number between 0 and 1, defines the sensitivity for aberration calling.')
	parser_test.add_argument('--blacklist',
							 type=str,
							 default=None,
							 help='Blacklist that masks regions in output, structure of header-less '
								  'file: chr...(/t)startpos(/t)endpos(/n)')
	parser_test.add_argument('--bed',
							 action='store_true',
							 help='Outputs tab-delimited .bed files, containing the most important information')
	parser_test.add_argument('--plot',
							 action='store_true',
							 help='Outputs .png plots')
	parser_test.set_defaults(func=tool_test)

	args = parser.parse_args(sys.argv[1:])

	logging.basicConfig(format='[%(levelname)s - %(asctime)s]: %(message)s',
						datefmt='%Y-%m-%d %H:%M:%S',
						level=getattr(logging, args.loglevel.upper(), None))
	logging.debug('args are: {}'.format(args))
	args.func(args)


if __name__ == '__main__':
	warnings.filterwarnings('ignore')
	main()