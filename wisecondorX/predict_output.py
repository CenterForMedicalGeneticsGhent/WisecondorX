# WisecondorX

import os

import numpy as np

from wisecondorX.overall_tools import exec_R, get_z_score, get_median_segment_variance, get_cpa
from pysam import VariantFile, VariantRecord

'''
Writes plots.
'''


def exec_write_plots(rem_input, results):
    json_plot_dir = os.path.abspath(rem_input['args'].outid + '_plot_tmp')
    json_dict = {
        'R_script': str('{}/include/plotter.R'.format(rem_input['wd'])),
        'ref_gender': str(rem_input['ref_gender']),
        'beta': str(rem_input['args'].beta),
        'zscore': str(rem_input['args'].zscore),
        'binsize': str(rem_input['binsize']),
        'n_reads': str(rem_input['n_reads']),
        'cairo': str(rem_input['args'].cairo),
        'results_r': results['results_r'],
        'results_w': results['results_w'],
        'results_c': results['results_c'],
        'ylim': str(rem_input['args'].ylim),
        'infile': str('{}.json'.format(json_plot_dir)),
        'out_dir': str('{}.plots'.format(rem_input['args'].outid)),
    }

    if rem_input["args"].add_plot_title:
        # Strip away paths from the outid if need be
        json_dict["plot_title"] = str(
            os.path.basename(rem_input["args"].outid))

    exec_R(json_dict)


'''
Calculates zz-scores, marks aberrations and
writes tables.
'''


def generate_output_tables(rem_input, results):
    if rem_input['args'].bed:
        _generate_bins_bed(rem_input, results)
    _generate_chr_statistics_file(rem_input, results)
    _generate_segments_and_aberrations(rem_input, results)

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
                r = 'nan'
            if z == 0:
                z = 'nan'
            feat_str = '{}:{}-{}'.format(chr_name, str(feat), str(feat + binsize - 1))
            row = [chr_name, feat, feat + binsize - 1, feat_str, r, z]
            bins_file.write('{}\n'.format('\t'.join([str(x) for x in row])))
            feat += binsize
    bins_file.close()


def _generate_segments_and_aberrations(rem_input, results):
    if rem_input['args'].bed:
        segments_file = open('{}_segments.bed'.format(rem_input['args'].outid), 'w')
        abberations_file = open('{}_aberrations.bed'.format(rem_input['args'].outid), 'w')
        segments_file.write('chr\tstart\tend\tratio\tzscore\n')
        abberations_file.write('chr\tstart\tend\tratio\tzscore\ttype\n')

    if rem_input['args'].vcf:
        vcf_file = VariantFile("{}.vcf.gz".format(rem_input['args'].outid), "w")
        prefix = _add_contigs(rem_input['args'].fai, vcf_file)
        _add_info(vcf_file)
        sample = rem_input['args'].sample if rem_input['args'].sample else rem_input['args'].outid.split("/")[-1]
        vcf_file.header.add_sample(sample)

    dup_count = 0
    del_count = 0
    cnv_count = 0

    for segment in results['results_c']:
        chr_name = str(segment[0] + 1)
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'
        start = int(segment[1] * rem_input['binsize'] + 1)
        stop = int(segment[2] * rem_input['binsize'])
        ratio = segment[4]
        zscore = segment[3]
        if rem_input['args'].bed: 
            row = [chr_name, start, stop, ratio, zscore]          
            segments_file.write('{}\n'.format('\t'.join([str(x) for x in row])))

        # output segments instead of abberations but add field that annotates which variant it is

        ploidy = 2
        if (chr_name == 'X' or chr_name == 'Y') and rem_input['ref_gender'] == 'M':
            ploidy = 1

        gain_or_loss = _define_gain_loss(segment, rem_input, ploidy)

        if rem_input['args'].vcf:
            if gain_or_loss == "gain":
                cnv_type = "DUP"
                dup_count += 1
                type_count = dup_count
            elif gain_or_loss == "loss":
                cnv_type = "DEL"
                del_count += 1
                type_count = del_count
            else:
                cnv_type = "CNV"
                cnv_count += 1
                type_count = cnv_count

            record: VariantRecord = vcf_file.new_record(
                contig=prefix + chr_name, 
                id="WisecondorX_{}_{}".format(cnv_type, type_count),
                start=start,
                stop=stop,
                alleles=('N', "<{}>".format(cnv_type))
            )
            record.info.update({
                'SVTYPE': 'CNV',
                'SVLEN': record.stop - record.start + 1,
                'ABB': False if not gain_or_loss else True,
                'SM': ratio,
                'ZS': zscore
            })
            record.samples[sample]['GT'] = (None,None) if ploidy == 2 else None
            vcf_file.write(record)

        if not gain_or_loss: continue
        if rem_input['args'].bed:
            abberations_file.write('{}\t{}\n'.format('\t'.join([str(x) for x in row]), gain_or_loss))

    if rem_input['args'].bed:
        segments_file.close()
        abberations_file.close()

    if rem_input['args'].vcf:
        vcf_file.close()

def _add_contigs(fai:str, vcf:VariantFile) -> str:
    """
    Adds the contigs to the VCF file
    """
    prefix = ""
    with open(fai, "r") as index:
        for line in index.readlines():
            if line.startswith("chr"): prefix = "chr"
            split_line = line.split("\t")
            vcf.header.add_line("##contig=<ID={},length={}>".format(split_line[0], split_line[1]))

    return prefix

def _add_info(vcf:VariantFile) -> None:
    """
    Adds the INFO fields to the VCF file
    """
    infos = [
        '##ALT=<ID=CNV,Description="Copy number variant region">',
        '##ALT=<ID=DEL,Description="Deletion relative to the reference">',
        '##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">',
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##INFO=<ID=SM,Number=1,Type=Float,Description="Linear copy ratio of the segment mean">',
        '##INFO=<ID=ZS,Number=1,Type=Float,Description="The z-score calculated for the current CNV">',
        '##INFO=<ID=ABB,Number=0,Type=Flag,Description="States that the CNV is an abberation">',
        '##FILTER=<ID=cnvQual,Description="CNV with quality below 10">',
        '##FILTER=<ID=cnvCopyRatio,Description="CNV with copy ratio within +/- 0.2 of 1.0">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    for info in infos:
        vcf.header.add_line(info)

def _define_gain_loss(segment, rem_input, ploidy = 2):
    """
    Define if the abberation is a gain or a loss
    """
    if rem_input['args'].beta is not None:
        abberation_cutoff = _get_aberration_cutoff(rem_input['args'].beta, ploidy)
        if float(segment[4]) > abberation_cutoff[1]:
            return "gain"
        elif float(segment[4]) < abberation_cutoff[0]:
            return "loss"
    elif isinstance(segment[3], str): return False
    else:
        if float(segment[3]) > rem_input['args'].zscore:
            return "gain"
        elif float(segment[3]) < - rem_input['args'].zscore:
            return "loss"
    return False

def _get_aberration_cutoff(beta, ploidy):
    loss_cutoff = np.log2((ploidy - (beta / 2)) / ploidy)
    gain_cutoff = np.log2((ploidy + (beta / 2)) / ploidy)
    return loss_cutoff, gain_cutoff


def _generate_chr_statistics_file(rem_input, results):
    stats_file = open('{}_statistics.txt'.format(rem_input['args'].outid), 'w')
    stats_file.write('chr\tratio.mean\tratio.median\tzscore\n')
    chr_ratio_means = [np.ma.average(results['results_r'][chr], weights=results['results_w'][chr])
                       for chr in range(len(results['results_r']))]
    chr_ratio_medians = [np.median([x for x in results['results_r'][chr] if x != 0])
                         for chr in range(len(results['results_r']))]

    results_c_chr = [[x, 0, rem_input['bins_per_chr'][x] - 1, chr_ratio_means[x]]
                     for x in range(len(results['results_r']))]

    msv = round(get_median_segment_variance(results['results_c'], results['results_r']), 5)
    cpa = round(get_cpa(results['results_c'], rem_input['binsize']), 5)
    chr_z_scores = get_z_score(results_c_chr, results)

    for chr in range(len(results['results_r'])):

        chr_name = str(chr + 1)
        if chr_name == '23':
            chr_name = 'X'
        if chr_name == '24':
            chr_name = 'Y'

        row = [chr_name,
               chr_ratio_means[chr],
               chr_ratio_medians[chr],
               chr_z_scores[chr]]

        stats_file.write('\t'.join([str(x) for x in row]) + '\n')

    stats_file.write('Gender based on --yfrac (or manually overridden by --gender): {}\n'
                     .format(str(rem_input['gender'])))

    stats_file.write('Number of reads: {}\n'
                     .format(str(rem_input['n_reads'])))

    stats_file.write('Standard deviation of the ratios per chromosome: {}\n'
                     .format(str(round(float(np.nanstd(chr_ratio_means)), 5))))

    stats_file.write('Median segment variance per bin (doi: 10.1093/nar/gky1263): {}\n'
                     .format(str(msv)))

    stats_file.write('Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4): {}\n'
                     .format(str(cpa)))

    stats_file.close()
