#!/usr/bin/env python

import argparse
from scipy.stats import norm
from wisecondorX.wisetools import *


def tool_convert(args):
    logging.info('Starting conversion')
    logging.info('Importing data ...')
    logging.info('Converting bam ... This might take a while ...')
    converted, qual_info = convert_bam(args.infile, binsize=args.binsize,
                                       min_shift=args.retdist, threshold=args.retthres)
    np.savez_compressed(args.outfile,
                        binsize=args.binsize,
                        sample=converted,
                        quality=qual_info)
    logging.info('Finished conversion')


def tool_newref(args):
    logging.info('Creating new reference')

    split_path = list(os.path.split(args.outfile))
    if split_path[-1][-4:] == '.npz':
        split_path[-1] = split_path[-1][:-4]
    base_path = os.path.join(split_path[0], split_path[1])

    # Add single thread information used for parallel processing functions
    args.prepfile = base_path + "_prep.npz"
    args.partfile = base_path + "_part"
    args.parts = args.cpus

    tool_newref_prep(args)

    # Use multiple cores if requested
    if args.cpus != 1:
        import concurrent.futures
        import copy
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpus) as executor:
            for part in xrange(1, args.parts + 1):
                if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
                    this_args = copy.copy(args)
                    this_args.part = [part, args.parts]
                    executor.submit(tool_newref_part, this_args)
            executor.shutdown(wait=True)
    else:
        for part in xrange(1, args.parts + 1):
            if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
                args.part = [part, args.parts]
                tool_newref_part(args)

    # Put it all together
    tool_newref_post(args)

    # Remove parallel processing temp data
    os.remove(args.prepfile)
    for part in xrange(1, args.parts + 1):
        os.remove(args.partfile + '_' + str(part) + '.npz')

    logging.info("Finished creating reference")


def tool_newref_prep(args):
    samples = []
    nreads = []
    binsizes = set()
    logging.info('Importing data ...')
    for infile in args.infiles:  # glob.glob(args.infiles):
        logging.info('Loading:{}'.format(infile))
        npzdata = np.load(infile)
        sample = npzdata['sample'].item()
        logging.info('binsize:{}'.format(int(npzdata['binsize'].item())))
        samples.append(scale_sample(sample, npzdata['binsize'].item(), args.binsize))
        nreads.append(sum([sum(sample[x]) for x in sample.keys()]))
        binsizes.add(npzdata['binsize'].item())

    if args.binsize is None and len(binsizes) != 1:
        logging.critical('There appears to be a mismatch in binsizes in your dataset: {} \n\t'
                         'Either remove the offending sample(s) or use a different --binsize'.format(binsizes))
        sys.exit()

    binsize = args.binsize
    if args.binsize is None:
        binsize = binsizes.pop()

    masked_data, chromosome_bins, mask = to_numpy_array(samples, args.gender)
    del samples
    masked_chrom_bins = [sum(mask[sum(chromosome_bins[:i]):sum(chromosome_bins[:i]) + x]) for i, x in
                         enumerate(chromosome_bins)]
    masked_chrom_bin_sums = [sum(masked_chrom_bins[:x + 1]) for x in range(len(masked_chrom_bins))]
    corrected_data, pca = train_pca(masked_data)
    np.savez_compressed(args.prepfile,
                        binsize=binsize,
                        chromosome_bins=chromosome_bins,
                        masked_data=masked_data,
                        mask=mask,
                        masked_chrom_bins=masked_chrom_bins,
                        masked_chrom_bin_sums=masked_chrom_bin_sums,
                        corrected_data=corrected_data,
                        pca_components=pca.components_,
                        pca_mean=pca.mean_)


def tool_newref_part(args):
    if args.part[0] > args.part[1]:
        logging.critical('Part should be smaller or equal to total parts:{} > {} is wrong'
                         .format(args.part[0], args.part[1]))
        sys.exit()
    if args.part[0] < 0:
        logging.critical('Part should be at least zero: {} < 0 is wrong'.format(args.part[0]))
        sys.exit()

    npzdata = np.load(args.prepfile)
    corrected_data = npzdata['corrected_data']
    masked_chrom_bins = npzdata['masked_chrom_bins']
    masked_chrom_bin_sums = npzdata['masked_chrom_bin_sums']

    logging.info('Creating reference ... This might take a while ...')
    indexes, distances = get_reference(corrected_data, masked_chrom_bins, masked_chrom_bin_sums, args.gender,
                                       select_ref_amount=args.refsize, part=args.part[0], split_parts=args.part[1])

    np.savez_compressed(args.partfile + '_' + str(args.part[0]) + '.npz',
                        indexes=indexes,
                        distances=distances)


def tool_newref_post(args):
    # Load prep file data
    npzdata = np.load(args.prepfile)
    masked_chrom_bins = npzdata['masked_chrom_bins']
    chromosome_bins = npzdata['chromosome_bins']
    mask = npzdata['mask']
    pca_components = npzdata['pca_components']
    pca_mean = npzdata['pca_mean']
    binsize = npzdata['binsize'].item()

    # Load and combine part file data
    big_indexes = []
    big_distances = []
    for part in xrange(1, args.parts + 1):  # glob.glob(args.infiles):
        infile = args.partfile + '_' + str(part) + '.npz'
        logging.info('Loading: {}'.format(infile))
        npzdata = np.load(infile)
        big_indexes.extend(npzdata['indexes'])
        big_distances.extend(npzdata['distances'])

        logging.info("{}, {}".format(part, npzdata['indexes'].shape))

    indexes = np.array(big_indexes)
    distances = np.array(big_distances)

    np.savez_compressed(args.outfile,
                        binsize=binsize,
                        indexes=indexes,
                        distances=distances,
                        chromosome_sizes=chromosome_bins,
                        mask=mask,
                        masked_sizes=masked_chrom_bins,
                        pca_components=pca_components,
                        pca_mean=pca_mean,
                        gender=args.gender)


def tool_test(args):
    logging.info('Starting CNA prediction')

    wc_dir = os.path.realpath(__file__)

    logging.info('Importing data ...')
    reference_file = np.load(args.reference)
    sample_file = np.load(args.infile)

    if not args.bed and not args.plot:
        logging.critical("No output format selected. \n\t"
                         "Select at least one of the supported output formats (--bed, --plot)")
        sys.exit()

    if args.beta <= 0 or args.beta > 1:
        logging.critical("Parameter beta should be a strictly positive number lower than 1")
        sys.exit()

    if args.alpha <= 0 or args.alpha > 1:
        logging.critical("Parameter alpha should be a strictly positive number lower than 1")
        sys.exit()

    if args.beta < 0.05:
        logging.warning("Parameter beta seems to be a bit low. \n\t"
                        "Have you read https://github.com/CenterForMedicalGeneticsGhent/wisecondorX#parameters on parameter optimization?")

    if args.alpha > 0.1:
        logging.warning("Parameter alpha seems to be a bit high. \n\t"
                        "Have you read https://github.com/CenterForMedicalGeneticsGhent/wisecondorX#parameters on parameter optimization?")

    # Reference data handling
    mask_list = []
    weights = get_weights(reference_file["distances"])
    gender = reference_file['gender']
    binsize = reference_file['binsize'].item()
    indexes = reference_file['indexes']
    distances = reference_file['distances']
    chromosome_sizes = reference_file['chromosome_sizes']
    mask = reference_file['mask']
    mask_list.append(mask)
    masked_sizes = reference_file['masked_sizes']
    masked_chrom_bin_sums = [sum(masked_sizes[:x + 1]) for x in range(len(masked_sizes))]

    pca_mean = reference_file['pca_mean']
    pca_components = reference_file['pca_components']

    del reference_file

    # Test sample data handling
    sample = sample_file['sample'].item()

    logging.info('Applying between-sample normalization ...')
    nreads = sum([sum(sample[x]) for x in sample.keys()])
    sample_bin_size = sample_file['binsize'].item()
    sample = scale_sample(sample, sample_bin_size, binsize)

    test_data = to_numpy_ref_format(sample, chromosome_sizes, mask, gender)
    test_data = apply_pca(test_data, pca_mean, pca_components)
    autosome_cutoff, allosome_cutoff = get_optimal_cutoff(distances, chromosome_sizes, args.maskrepeats)

    z_threshold = norm.ppf(0.975)  # two-tailed test

    logging.info('Applying within-sample normalization ...')
    test_copy = np.copy(test_data)
    results_z, results_r, ref_sizes, std_dev_avg = repeat_test(test_copy, indexes, distances,
                                                               masked_sizes, masked_chrom_bin_sums,
                                                               autosome_cutoff, allosome_cutoff, z_threshold, 5)
    # Get rid of infinite values caused by having no reference bins or only zeros in the reference
    infinite_mask = (ref_sizes >= args.minrefbins)
    mask_list.append(infinite_mask)
    cleaned_r = results_r[infinite_mask]
    cleaned_z = results_z[infinite_mask]
    cleaned_weights = weights[infinite_mask]
    cleaned_bin_sums = [np_sum(infinite_mask[:val]) for val in masked_chrom_bin_sums]
    cleaned_bins = [cleaned_bin_sums[i] - cleaned_bin_sums[i - 1] for i in range(1, len(cleaned_bin_sums))]
    cleaned_bins.insert(0, cleaned_bin_sums[0])


    results_z = []
    results_r = []
    results_w = []
    inflated_z = inflate_array_multi(cleaned_z, mask_list)
    inflated_r = inflate_array_multi(cleaned_r, mask_list)
    inflated_w = inflate_array_multi(cleaned_weights, mask_list)
    for chrom in xrange(len(chromosome_sizes)):
        chrom_data = inflated_z[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
        results_z.append(chrom_data)
        chrom_data = inflated_r[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
        results_r.append(chrom_data)
        chrom_data = inflated_w[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
        results_w.append(chrom_data)

    # log2
    for chrom in xrange(len(chromosome_sizes)):
        results_r[chrom] = np.log2(results_r[chrom])

    # Apply blacklist
    if args.blacklist:
        apply_blacklist(args, binsize, results_r, results_z, results_w, sample, gender)

    # Make R interpretable
    results_r = [x.tolist() for x in results_r]
    results_z = [x.tolist() for x in results_z]
    results_w = [x.tolist() for x in results_w]
    nchrs = len(results_r)
    for c in range(nchrs):
        for i, rR in enumerate(results_r[c]):
            if not np.isfinite(rR):
                results_r[c][i] = 0
                results_z[c][i] = 0
                results_w[c][i] = 0

    logging.info('Obtaining CBS segments ...')
    cbs_calls = cbs(args, results_r, results_z, results_w, gender, wc_dir)

    json_out = {'binsize': binsize,
                'results_r': results_r,
                'results_z': results_z,
                'results_w' : results_w,
                'threshold_z': z_threshold.tolist(),
                'asdef': std_dev_avg.tolist(),
                'aasdef': (std_dev_avg * z_threshold).tolist(),
                'nreads': nreads,
                'cbs_calls': cbs_calls}

    # Save txt: optional
    if args.bed:
        logging.info('Writing tables ...')
        generate_txt_output(args, json_out)

    # Create plots: optional
    if args.plot:
        logging.info('Writing plots ...')
        write_plots(args, json_out, wc_dir, gender)

    logging.info("Finished prediction")


def get_gender(args):
    npzfile = np.load(args.infile)
    sample = npzfile['sample'].item()
    non_y = float(sum([np.sum(sample[str(chr)]) for chr in range(1, 24)]))
    y = float(np.sum(sample["24"]))
    permille_y = y / (non_y + y) * 1000.0
    if permille_y > args.cutoff:
        print("male")
    else:
        print("female")


def main():
    parser = argparse.ArgumentParser(description="wisecondorX")
    parser.add_argument('--loglevel',
                        type=str,
                        default='INFO',
                        choices=['info', 'warning', 'debug', 'error', 'critical'])
    subparsers = parser.add_subparsers()

    # File conversion
    parser_convert = subparsers.add_parser('convert',
                                           description='Convert and filter a .bam file to a .npz',
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_convert.add_argument('infile',
                                type=str,
                                help='Bam input file for conversion')
    parser_convert.add_argument('outfile',
                                type=str,
                                help='Output npz file')
    parser_convert.add_argument('--binsize',
                                type=float,
                                default=5e3,
                                help='Bin size (bp).')
    parser_convert.add_argument('--retdist',
                                type=int,
                                default=4,
                                help='Maximum amount of base pairs difference between sequential reads '
                                     'to consider them part of the same tower.')
    parser_convert.add_argument('--retthres',
                                type=int,
                                default=4,
                                help='Threshold for when a group of reads is considered a tower and will be removed.')
    parser_convert.set_defaults(func=tool_convert)

    # Find gender
    parser_gender = subparsers.add_parser('gender',
                                          description='Predict the gender of a sample',
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gender.add_argument('infile',
                               type=str,
                               help='.npz input file')
    parser_gender.add_argument('--cutoff',
                               type=float,
                               default=3.5,
                               help='Y-read permille cut-off. Below is female, above is male.')
    parser_gender.set_defaults(func=get_gender)

    # New reference creation
    parser_newref = subparsers.add_parser('newref',
                                          description='Create a new reference using healthy reference samples',
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_newref.add_argument('infiles',
                               type=str,
                               nargs='+',
                               help='Path to all reference data files (e.g. path/to/reference/*.npz)')
    parser_newref.add_argument('outfile',
                               type=str,
                               help='Path and filename for the reference output (i.e. ./reference/myref.npz)')
    parser_newref.add_argument('--refsize',
                               type=int,
                               default=300,
                               help='Amount of reference locations per target.')
    parser_newref.add_argument('--binsize',
                               type=int,
                               default=1e5,
                               help='Scale samples to this binsize, multiples of existing binsize only.')
    parser_newref.add_argument('--gender',
                               type=str,
                               default="F",
                               choices=["F", "M"],
                               help='Gender of the reference .npz input files.')
    parser_newref.add_argument('--cpus',
                               type=int,
                               default=1,
                               help='Use multiple cores to find reference bins')
    parser_newref.set_defaults(func=tool_newref)

    # Find CNAs
    parser_test = subparsers.add_parser('predict',
                                        description='Find copy number aberrations',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_test.add_argument('infile',
                             type=str,
                             help='Sample.npz of which the CNAs will be predicted')
    parser_test.add_argument('reference',
                             type=str,
                             help='Reference as previously created')
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
                             help='P-value cut-off for calling a CBS breakpoint.')
    parser_test.add_argument('--beta',
                             type=float,
                             default=0.1,
                             help='Number between 0 and 1, defines the sensitivity for aberration calling.')
    parser_test.add_argument('--blacklist',
                             type=str,
                             default=None,
                             help='Blacklist that masks regions in output, structure of header-less '
                                  'file: chrX(/t)startpos(/t)endpos(/n)')
    parser_test.add_argument('--bed',
                             action="store_true",
                             help='Outputs tab-delimited .bed files, containing the most important information')
    parser_test.add_argument('--plot',
                             action="store_true",
                             help='Outputs .png plots')
    parser_test.set_defaults(func=tool_test)

    args = parser.parse_args(sys.argv[1:])

    # configure logging
    logging.basicConfig(format='[%(levelname)s - %(asctime)s]: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=getattr(logging, args.loglevel.upper(), None))
    logging.debug("args are: {}".format(args))
    args.func(args)


if __name__ == '__main__':
    import warnings

    warnings.filterwarnings("ignore")
    main()
