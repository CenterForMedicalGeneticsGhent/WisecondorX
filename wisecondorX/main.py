#!/usr/bin/env python

import argparse
from scipy.stats import norm
from wisecondorX.wisetools import *


def tool_convert(args):
    logging.info('Starting conversion')
    logging.info('Importing data ...')
    logging.info('Converting bam ... This might take a while ...')
    sample, qual_info = convert_bam(args.infile, binsize=args.binsize,
                                       min_shift=args.retdist, threshold=args.retthres)
    if args.gender:
        gender = args.gender
        logging.info('Gender {}'.format(gender))
    else:
        logging.info('Predicting gender ...')
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

    # Add single thread information used for parallel processing functions
    args.basepath = base_path
    args.prepfile = base_path + "_prep.npz"
    args.partfile = base_path + "_part"

    samples = []
    genders = []
    binsizes = set()
    logging.info('Importing data ...')
    for infile in args.infiles:
        logging.info('Loading: {}'.format(infile))
        npzdata = np.load(infile)
        sample = npzdata['sample'].item()
        binsize = npzdata['binsize'].item()
        gender = npzdata['gender'].item()
        logging.info('Binsize: {} | Gender : {}'.format(int(binsize), gender))

        scaled_sample = scale_sample(sample, binsize, args.binsize)
        corrected_sample = gender_correct(scaled_sample, gender)

        genders.append(gender)
        samples.append(corrected_sample)
        binsizes.add(binsize)

    if args.binsize is None and len(binsizes) != 1:
        logging.critical('There appears to be a mismatch in binsizes in your dataset: {} \n\t'
                         'Either remove the offending sample(s) or use a different --binsize'.format(binsizes))
        sys.exit()

    samples = np.array(samples)

    outfiles = []
    if len(genders) > 9:
        logging.info('Starting autosomal reference creation ...')
        args.tmpoutfile = args.basepath + ".tmp.A.npz"
        outfiles.append(args.tmpoutfile)
        tool_newref_prep(args, samples, np.array(genders), "A")
        logging.info('This might take a while ...')
        tool_newref_main(args, args.cpus)
    else:
        logging.critical('Provide at least 10 samples to enable the generation of a reference.')
        sys.exit()

    if genders.count("F") > 4:
        logging.info('Starting female gonosomal reference creation ...')
        args.tmpoutfile = args.basepath + ".tmp.F.npz"
        outfiles.append(args.tmpoutfile)
        tool_newref_prep(args, samples, np.array(genders), "F")
        logging.info('This might take a while ...')
        tool_newref_main(args, 1)
    else:
        logging.warning('Provide at least 5 female samples to enable normalization of female gonosomes.')

    if genders.count("M") > 4:
        logging.info('Starting male gonosomal reference creation ...')
        args.tmpoutfile = args.basepath + ".tmp.M.npz"
        outfiles.append(args.tmpoutfile)
        tool_newref_prep(args, samples, np.array(genders), "M")
        tool_newref_main(args, 1)
    else:
        logging.warning('Provide at least 5 male samples to enable normalization of male gonosomes. '
                        'If these are of no interest (e.g. NIPT), ignore this warning.')

    tool_newref_merge(args, outfiles)
    logging.info("Finished creating reference")


def tool_newref_main(args, cpus):

    # Use multiple cores if requested
    if cpus != 1:
        import concurrent.futures
        import copy
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpus) as executor:
            for part in xrange(1, cpus + 1):
                if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
                    this_args = copy.copy(args)
                    this_args.part = [part, cpus]
                    executor.submit(tool_newref_part, this_args)
            executor.shutdown(wait=True)
    else:
        for part in xrange(1, cpus + 1):
            if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
                args.part = [part, cpus]
                tool_newref_part(args)

    # Put it together
    tool_newref_post(args, cpus)

    # Remove parallel processing temp data
    os.remove(args.prepfile)
    for part in xrange(1, cpus + 1):
        os.remove(args.partfile + '_' + str(part) + '.npz')


def tool_newref_prep(args, samples, genders, gender):

    if gender == "A":
        masked_data, chromosome_bins, mask = to_numpy_array(samples, range(1, 23))
    elif gender == "F":
        masked_data, chromosome_bins, mask = to_numpy_array(samples[genders == gender], range(1, 24))
    else:
        masked_data, chromosome_bins, mask = to_numpy_array(samples[genders == gender], range(1, 25))

    del samples
    masked_chrom_bins = [sum(mask[sum(chromosome_bins[:i]):sum(chromosome_bins[:i]) + x]) for i, x in
                         enumerate(chromosome_bins)]
    masked_chrom_bin_sums = [sum(masked_chrom_bins[:x + 1]) for x in range(len(masked_chrom_bins))]
    corrected_data, pca = train_pca(masked_data)
    np.savez_compressed(args.prepfile,
                        binsize=args.binsize,
                        chromosome_bins=chromosome_bins,
                        masked_data=masked_data,
                        mask=mask,
                        masked_chrom_bins=masked_chrom_bins,
                        masked_chrom_bin_sums=masked_chrom_bin_sums,
                        corrected_data=corrected_data,
                        gender=gender,
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

    indexes, distances = get_reference(corrected_data, masked_chrom_bins, masked_chrom_bin_sums,
                                       select_ref_amount=args.refsize, part=args.part[0], split_parts=args.part[1])

    np.savez_compressed(args.partfile + '_' + str(args.part[0]) + '.npz',
                        indexes=indexes,
                        distances=distances)


def tool_newref_post(args, cpus):
    # Load prep file data
    npzdata = np.load(args.prepfile)
    masked_chrom_bins = npzdata['masked_chrom_bins']
    chromosome_bins = npzdata['chromosome_bins']
    mask = npzdata['mask']
    pca_components = npzdata['pca_components']
    pca_mean = npzdata['pca_mean']
    binsize = npzdata['binsize'].item()
    gender = npzdata['gender'].item()

    # Load and combine part file data
    big_indexes = []
    big_distances = []
    for part in xrange(1, cpus + 1):  # glob.glob(args.infiles):
        infile = args.partfile + '_' + str(part) + '.npz'
        npzdata = np.load(infile)
        big_indexes.extend(npzdata['indexes'])
        big_distances.extend(npzdata['distances'])

    indexes = np.array(big_indexes)
    distances = np.array(big_distances)

    np.savez_compressed(args.tmpoutfile,
                        binsize=binsize,
                        indexes=indexes,
                        distances=distances,
                        chromosome_sizes=chromosome_bins,
                        mask=mask,
                        masked_sizes=masked_chrom_bins,
                        pca_components=pca_components,
                        pca_mean=pca_mean,
                        gender=gender)


def tool_newref_merge(args, outfiles):
    final_ref = {"has_female" : False, "has_male" : False}
    for file_id in outfiles:
        npz_file = np.load(file_id)
        gender = str(npz_file['gender'])
        for component in [x for x in npz_file.keys() if x != "gender"]:
            if gender == "F":
                final_ref["has_female"] = True
                final_ref[str(component) + ".F"] = npz_file[component]
            elif gender == "M":
                final_ref["has_male"] = True
                final_ref[str(component) + ".M"] = npz_file[component]
            else:
                final_ref[str(component)] = npz_file[component]
        os.remove(file_id)
    np.savez_compressed(args.outfile, **final_ref)


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
    binsize = reference_file['binsize'].item()
    indexes = reference_file['indexes']
    distances = reference_file['distances']
    chromosome_sizes = reference_file['chromosome_sizes']
    mask = reference_file['mask']
    masked_sizes = reference_file['masked_sizes']
    masked_chrom_bin_sums = [sum(masked_sizes[:x + 1]) for x in range(len(masked_sizes))]
    pca_mean = reference_file['pca_mean']
    pca_components = reference_file['pca_components']
    ref_has_male = reference_file['has_male']
    ref_has_female = reference_file['has_female']

    # Test sample data handling
    sample = sample_file['sample'].item()
    nreads = sum([sum(sample[x]) for x in sample.keys()])
    actual_gender = sample_file['gender']
    reference_gender = str(actual_gender)
    sample = gender_correct(sample, actual_gender)

    logging.info('Applying between-sample normalization autosomes...')
    sample_bin_size = sample_file['binsize'].item()
    sample = scale_sample(sample, sample_bin_size, binsize)

    test_data = to_numpy_ref_format(sample, chromosome_sizes, mask)
    test_data = apply_pca(test_data, pca_mean, pca_components)
    cutoff = get_optimal_cutoff(distances, args.maskrepeats, 0)

    z_threshold = norm.ppf(0.975)  # two-tailed test

    logging.info('Applying within-sample normalization autosomes...')
    test_copy = np.copy(test_data)
    results_z, results_r, ref_sizes = repeat_test(test_copy, indexes, distances,
                                                masked_sizes, masked_chrom_bin_sums,
                                                cutoff, z_threshold, 5)

    if not ref_has_male and actual_gender == "M":
        logging.warning('This sample is male, whilst the reference is created with fewer than 5 males. '
                        'The female gonosomal reference will be used for X predictions. Note that these might '
                        'not be accurate. If the latter is desired, create a new reference and include more '
                        'male samples.')
        reference_gender = "F"

    elif not ref_has_female and actual_gender == "F":
        logging.warning('This sample is female, whilst the reference is created with fewer than 5 females. '
                        'The male gonosomal reference will be used for XY predictions. Note that these might '
                        'not be accurate. If the latter is desired, create a new reference and include more '
                        'female samples.')
        reference_gender = "M"

    results_z, results_r, ref_sizes, weights, chromosome_sizes, mask, masked_sizes, masked_chrom_bin_sums = \
        append_objects_with_gonosomes(args, reference_gender, sample, reference_file,
                                  z_threshold,results_z, results_r,
                                  ref_sizes, weights, mask, masked_sizes)
    del reference_file
    mask_list.append(mask)

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
        apply_blacklist(args, binsize, results_r, results_z, results_w)

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
    cbs_calls = cbs(args, results_r, results_z, results_w, reference_gender, wc_dir)

    out_dict = {'binsize': binsize,
                'results_r': results_r,
                'results_z': results_z,
                'results_w' : results_w,
                'nreads': nreads,
                'cbs_calls': cbs_calls,
                'actual_gender': str(actual_gender),
                'reference_gender': str(reference_gender),
                'beta': str(args.beta)}

    # Save txt: optional
    if args.bed:
        logging.info('Writing tables ...')
        generate_txt_output(args, out_dict)

    # Create plots: optional
    if args.plot:
        logging.info('Writing plots ...')
        write_plots(args, out_dict, wc_dir)

    logging.info("Finished prediction")


def output_gender(args):
    npzdata = np.load(args.infile)
    gender = str(npzdata['gender'])
    if gender == "M":
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
                                choices=["F", "M"],
                                help='Gender of the case. If not given, WisecondorX will predict it')
    parser_convert.add_argument('--gonmapr',
                                type=float,
                                default=2,
                                help='The gonosomal mappabality ratio between X and Y. Concerning short single-end '
                                     'read mapping, a Y bin is two times (default) less mappable compared to an X bin. '
                                     'Used to predict gender')
    parser_convert.set_defaults(func=tool_convert)

    # Get gender
    parser_gender = subparsers.add_parser('gender',
                                           description='Returns the gender of a .npz resulting from convert',
                                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_gender.add_argument('infile',
                                type=str,
                                help='.npz input file')
    parser_gender.set_defaults(func=output_gender)

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

    # Find CNAs
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