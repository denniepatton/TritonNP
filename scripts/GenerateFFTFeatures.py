#!/usr/bin/env python
# Robert Patton, rpatton@fredhutch.org (Ha Lab)
# v4.1, 03/08/2022 : cleaned up for GitHub

# Takes as input a bam file, associated GC bias, and a bed file with regions of interest, and returns a feature
# matrix with two features: "NPS" (Nucleosome Phasing Score) and "mean-period"
# a range of fragment sizes to be considered can also be specified

# Method: pileup of raw coverage is run through a fast Fourier transform and high-frequency signal corresponding
# to periodic signal less than 146 bp (smaller than the minimum size of nucleosome-originating periodicity) is
# filtered out. Peaks are then called on the reconstructed signal. mean-period is the number of peaks divided by
# the region of interest length, and phased-component is the ratio of typical to compact signal component amplitudes

# illustrative figures can be output by specifying them in the "plot_list" file


import os
import sys
import pysam
import random
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.fft import rfft, rfftfreq, irfft
from multiprocessing import Pool
from functools import partial


def get_gc_bias_dict(bias_path):
    if bias_path is not None:
        bias = pd.read_csv(bias_path, sep='\t')
        bias['smoothed_GC_bias'] = np.where(bias['smoothed_GC_bias'] < 0.05, np.nan, bias['smoothed_GC_bias'])
        bias = bias[['length', 'num_GC', 'smoothed_GC_bias']]
        bias = bias.set_index(['num_GC', 'length']).unstack()
        bias = bias.to_dict()
        bias2 = {}
        for key in bias.keys():
            length = key[1]
            bias2[length] = {}
            for num_GC in range(0, length + 1):
                temp_bias = bias[key][num_GC]
                bias2[length][num_GC] = temp_bias
        bias = bias2
        del bias2
        return bias
    else:
        return None


def local_peaks(ys, xs):
    max_values, max_indices = [], []
    ysl = len(ys) - 1
    for i, y in enumerate(ys):
        if 0 < i < ysl and ys[i - 1] <= y and y > ys[i + 1] and xs[i] > 0:  # local maxima
            max_values.append(y)
            max_indices.append(i)
    return max_values, max_indices


def fft_info(bed_region, params):
    bam_path, sample, out_direct, frag_range, gc_bias, ref_seq_path, to_plot, map_q = params
    bam = pysam.AlignmentFile(bam_path, 'rb')
    ref_seq = pysam.FastaFile(ref_seq_path)
    # theoretical nucleosome period is 146 bp -> f = 0.0068493
    freq_max = 0.0068493
    # range_1: T = 150-180 bp (f = 0.00556 - 0.00667)
    # range_2: T = 180-210 bp (f = 0.00476 - 0.00555)
    low_1, high_1 = 0.00556, 0.00667
    low_2, high_2 = 0.00476, 0.00555
    bed_tokens = bed_region.strip().split('\t')
    # below I am adding a 500 bp buffer to better fit the window itself
    start_pos = int(bed_tokens[1]) - 500
    stop_pos = int(bed_tokens[2]) + 500
    site = str(bed_tokens[3])
    depth = [0] * (stop_pos - start_pos)
    roi_length = len(depth)
    # process all fragments falling inside the ROI
    segment_reads = bam.fetch(bed_tokens[0], start_pos, stop_pos)
    for read in segment_reads:
        fragment_length = read.template_length
        if frag_range[0] <= np.abs(fragment_length) <= frag_range[1] and read.is_paired and read.\
                mapping_quality >= map_q and not read.is_duplicate and not read.is_qcfail:
            read_start = read.reference_start - start_pos
            if read.is_reverse and fragment_length < 0:
                read_length = read.reference_length
                fragment_start = read_start + read_length + fragment_length
                fragment_end = read_start + read_length
            elif not read.is_reverse and fragment_length > 0:
                fragment_start = read_start
                fragment_end = read_start + fragment_length
            else:
                continue
            # now do GC bias correction:
            if gc_bias is not None:
                fragment_seq = (ref_seq.fetch(bed_tokens[0], fragment_start + start_pos, fragment_end + start_pos)).\
                    upper()
                fragment_seq = list(fragment_seq.replace('T', '0').replace('A', '0').replace('C', '1').
                                    replace('G', '1').replace('N', str(np.random.randint(0, 2))))
                fragment_gc_content = sum([int(m) for m in fragment_seq])
                fragment_bias = gc_bias[np.abs(fragment_length)][fragment_gc_content]
            else:
                fragment_bias = 1
            fragment_cov = np.array(range(fragment_start, fragment_end + 1))
            for index in np.clip(fragment_cov, 0, roi_length - 1):
                if 0.05 < fragment_bias < 10:
                    depth[index] += 1 / fragment_bias
    mean_depth = np.mean(depth[500:-500])
    fourier = rfft(np.asarray(depth))
    freqs = rfftfreq(roi_length)
    range_1 = [idx for idx, val in enumerate(freqs) if low_1 <= val <= high_1]
    range_2 = [idx for idx, val in enumerate(freqs) if low_2 <= val <= high_2]
    primary_amp_1 = round(np.mean(np.abs(fourier[range_1])) / roi_length, 4)
    primary_amp_2 = round(np.mean(np.abs(fourier[range_2])) / roi_length, 4)
    if primary_amp_1 > 0 and primary_amp_2 > 0:
        amp_ratio = primary_amp_2 / primary_amp_1
    else:
        amp_ratio = np.nan
    test_freqs = [idx for idx, val in enumerate(freqs) if 0 < val <= freq_max]  # frequencies in filter
    clean_fft = [0] * len(fourier)
    try:
        clean_fft[test_freqs[0]:test_freqs[-1]] = fourier[test_freqs[0]:test_freqs[-1]]
    except IndexError:  # not enough data to construct a reasonable Fourier
        return site, np.nan, np.nan, np.nan
    inverse_signal = irfft(clean_fft)  # reconstruct signal
    primary_signal = inverse_signal + len(inverse_signal) * [mean_depth]  # add in base component
    # Remove the buffer regions:
    depth = depth[500:-500]
    primary_signal = primary_signal[500:-500]
    max_values, peaks = local_peaks(primary_signal, depth)
    if len(max_values) < 1:
        nuc_period = np.nan
    else:
        nuc_period = round(roi_length / len(peaks), 2)
    if site in to_plot and mean_depth > 1:  # plot selected regions
        clean_1 = [0] * len(fourier)
        clean_2 = [0] * len(fourier)
        clean_1[range_1[0]:range_1[-1]] = fourier[range_1[0]:range_1[-1]]
        clean_2[range_2[0]:range_2[-1]] = fourier[range_2[0]:range_2[-1]]
        inverse_1 = irfft(clean_1)[500:-500]
        inverse_2 = irfft(clean_2)[500:-500]
        fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 8))
        fig.suptitle(bed_tokens[3] + ' Depth and FFT', fontsize=20)
        ax1.fill_between(list(range(len(depth))), 0, depth, label='Coverage', color='gray')
        ax1.set_xlim(0, len(depth))
        ax1.set_ylim(0, 3 * mean_depth)
        ax1.set_ylabel('Depth', fontsize=20)
        ax1.set_xlabel('Gene Coordinates', fontsize=20)
        if len(peaks) > 1:
            ax1.axvline(x=peaks[0], color='gray', label='Local Peaks', linestyle='dashed', alpha=0.2)
            for loc in peaks[1:]:
                ax1.axvline(x=loc, color='gray', linestyle='dashed', alpha=0.2)
        ax1.plot(inverse_1 + len(inverse_1) * [mean_depth * 0.5], color='#C03830',
                 label='Nucleosome Component (Avg Amp = ' + str(primary_amp_1) + ')')
        ax1.plot(inverse_2 + len(inverse_2) * [mean_depth * 0.5], color='#49B022',
                 label='Nucleosome + Linker Component (Avg Amp = ' + str(primary_amp_2) + ')')
        ax1.plot(primary_signal, color='#8927D6FA', label='Phased Component')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=20)
        ax2.plot(freqs, np.abs(fourier) / roi_length)
        ax2.axvspan(freqs[test_freqs[0]], freqs[test_freqs[-1]], facecolor='#8927D6', alpha=0.2)
        ax2.axvline(x=(1/nuc_period), color='gray', linestyle='dashed',
                    label='Peak-Based f: ' + str(round(1/nuc_period, 6)) + ' /bp (T = ' + str(nuc_period) + ' bp)')
        ax2.set_xlim(0, 0.01)
        ax2.set_ylim(0, 10)
        ax2.set_xlabel('Frequency (/bp)', fontsize=20)
        ax2.set_ylabel('Amplitude', fontsize=20)
        ax2.legend(loc=2, fontsize=20)
        fig.savefig(out_direct + '/' + sample + '_' + site + '.pdf', bbox_inches="tight")
        plt.close()
    return site, nuc_period, amp_ratio, mean_depth


def main():
    # parse command line arguments:
    parser = argparse.ArgumentParser(description='\n### GenerateFragmentFeatures.py ### Combine a bam coverage file and'
                                                 ' a bed file with regions of interest to produce a feature matrix '
                                                 'with samples as rows, and features as columns. Features are '
                                                 'NPS and mean-period. See code for details and/or to '
                                                 'change which regions get plotted.')
    parser.add_argument('-n', '--sample_name', help='sample identifier', required=True)
    parser.add_argument('-i', '--input', help='bam file', required=True)
    parser.add_argument('-b', '--bias', help='GC bias file (from Griffin)', default=None)
    parser.add_argument('-a', '--annotation', help='bed file with regions of interest', required=True)
    parser.add_argument('-g', '--reference_genome', help='path to reference genome', required=True)
    parser.add_argument('-r', '--results_dir', help='directory for results', required=True)
    parser.add_argument('-p', '--plot_list', help='list of genes/regions to generate plots for', required=True)
    parser.add_argument('-q', '--map_quality', help='minimum mapping quality', type=int, default=20)
    parser.add_argument('-f', '--size_range', help='fragment size range to use (bp)', nargs=2, type=int, default=(15, 500))
    parser.add_argument('-c', '--cpus', help='cpu available for parallel processing', type=int, required=True)
    args = parser.parse_args()

    print('Loading input files . . .')

    sample_name = args.sample_name
    bam_path = args.input
    bias_path = args.bias
    bed_path = args.annotation
    ref_seq_path = args.reference_genome
    results_dir = args.results_dir
    plot_list = args.plot_list
    map_q = args.map_quality
    size_range = args.size_range
    cpus = args.cpus

    print('\narguments provided:')
    print('\tsample_name = "' + sample_name + '"')
    print('\tbam_path = "' + bam_path + '"')
    print('\tGC_bias_path = "' + bias_path + '"')
    print('\tref_seq_path = "' + ref_seq_path + '"')
    print('\tsites_bed = "' + bed_path + '"')
    print('\tplot_list = "' + plot_list + '"')
    print('\tresults_dir = "' + results_dir + '"')
    print('\tsize_range =', size_range)
    print('\tmap_q =', map_q)
    print('\tCPUs =', cpus)
    print('\n')
    sys.stdout.flush()

    gc_bias = get_gc_bias_dict(bias_path)
    to_plot = pd.read_table(plot_list, header=None)[0].tolist()
    sites = [region for region in open(bed_path, 'r')]
    random.shuffle(sites)
    params = [bam_path, sample_name, results_dir, size_range, gc_bias, ref_seq_path, to_plot, map_q]

    print('Running fft_info on ' + str(len(sites)) + ' regions . . .')

    with Pool(cpus) as pool:
        results = list(pool.imap_unordered(partial(fft_info, params=params), sites, len(sites) // cpus))

    print('Merging results . . .')

    fm = {sample_name: {'Sample': sample_name}}
    for result in results:
        fm[sample_name][result[0] + '_peak-based-period'] = result[1]
        fm[sample_name][result[0] + '_amplitude-ratio'] = result[2]
        fm[sample_name][result[0] + '_mean-depth'] = result[3]
    df = pd.DataFrame(fm).transpose()

    out_file = results_dir + '/' + sample_name + '_PhasingFM.tsv'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    df.to_csv(out_file, sep='\t')

    print('Finished.')


if __name__ == "__main__":
    main()
