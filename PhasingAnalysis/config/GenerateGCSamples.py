#!/usr/bin/python
# Robert Patton, rpatton@fredhutch.org

# simple command line tool for making the sample files used by Griffin and GenerateFFTFeatures

# v3.0, 11/02/2021

import getopt
import sys
import os
from glob import glob
from collections import defaultdict


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:b:v", ["help", "output=", "bams=", "biases=", "verbose"])
    except getopt.GetoptError as err:
        print(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    if not opts:
        usage()
        sys.exit(2)
    output = "samples.yaml"  # default name
    verbose, bias_direct, bam_direct = False, None, None
    for o, a in opts:
        if o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-i", "--bams"):
            bam_direct = a
        elif o in ("-b", "--biases"):
            bias_direct = a
        else:
            assert False, "unhandled option"

    bam_paths = [y for x in os.walk(bam_direct) for y in glob(os.path.join(x[0], '*.bam'))]
    bam_dict = dict(zip(map(os.path.basename, bam_paths), bam_paths))
    bam_dict = {x.replace('.bam', ''): v for x, v in bam_dict.items()}

    bias_paths = [y for x in os.walk(bias_direct) for y in glob(os.path.join(x[0], '*.GC_bias.txt'))]
    bias_dict = dict(zip(map(os.path.basename, bias_paths), bias_paths))
    bias_dict = {x.replace('.GC_bias.txt', ''): v for x, v in bias_dict.items()}

    sample_dict = defaultdict(list)
    for d in (bam_dict, bias_dict):
        for key, value in d.items():
            sample_dict[key].append(value)

    total = 0
    if verbose:
        print("\nPrinting the following samples to " + output + " . . .\n")
    with open(output, "w") as out:
        out.write("samples:\n")
        for sample, info in sample_dict.items():
            if len(info) == 2:  # vs 3, before
                out.write("  " + sample + ":\n")
                out.write("    bam: " + info[0] + "\n")
                out.write("    GC_bias: " + info[1] + "\n")
                if verbose:
                    print(sample)
                    total += 1
            elif verbose:
                print("### WARNING ### " + sample + " not written: missing inputs\n"
                      "### provided inputs:\n")
                for item in info:
                    print("### " + item + "\n")
    if verbose:
        print("\nFinished. " + str(total) + " samples written.\n")


def usage():
    print("\n### GenerateGCSamples.py ###\n"
          "Combines bam and bias information from multiple directories to generate\n"
          "a samples (.yaml) file for use with the Griffin-like pipeline. bam and bias files come from GC correction.\n"
          "All files for a single sample should contain a common naming convention, and this script will search input\n"
          " folders recursively by default.\n"
          "\n"
          "Usage: GenerateGCSamples.py [options] -i|--bams <bams_direct> -b|--biases <biases_direct>\n"
          "Options:\n"
          "  -h|--help          show this usage message\n"
          "  -o|--output        specify output file (default: samples.yaml)\n"
          "  -v|--verbose       run in verbose mode: print descriptive dialogue\n")


if __name__ == "__main__":
    main()
