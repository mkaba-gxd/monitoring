import os
import sys
import re
import glob
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from cyvcf2 import VCF
from .func import *

use_chrom = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
use_columns = ["CHROM", "POS", "REF", "ALT"]

def load_vcf(vcf_path):

    with open(vcf_path) as f:
        skip_rows = sum(1 for line in f if line.startswith('#'))

    data = pd.read_csv(vcf_path, sep='\t', skiprows=skip_rows, header=None, index_col=False, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','ATTR'])
    data = data[use_columns].copy()
    return data

def search_vcf(vcf_path, chrom, position, window=0):

    data = load_vcf(vcf_path)
    data = data[ data['CHROM']==chrom ]
    if window == 0:
        data = data[data['POS']==int(position)].reset_index(drop=True)
    else :
        data = data[(int(position)-window<=data['POS']) & (data['POS']<=int(position)+window)].reset_index(drop=True)

    return data

def run_snv(args):

    sample = args.sample
    locus = args.position.split(':')
    window = args.window
    anal_dir = args.directory

    if len(locus) != 2:
        init('input value error: --position')
    elif not locus[0] in use_chrom :
        init('input value error: No chromosome ' + locus[0])
    elif not locus[1].isnumeric() :
        init('The position should be entered as an integer.')

    subDir = batch(sample, anal_dir, 'eWES')
    if subDir is None : init("No registration in database")

    # mutect2
    file = os.path.join(anal_dir, 'eWES', subDir, sample, 'SNV', 'somatic', 'mutect2', sample + '.mutect2.cleaned.vcf')
    if not os.path.isfile(file) :
        print("mutect2: Intermediate file does not exist.")
    else :
        data = search_vcf(file, locus[0], locus[1], window)
        if data.shape[0] == 0 :
            print('mutect2: No matching data found.')
        else :
            print('[mutect2]')
            print(data.to_string(index=False))

    # lofreq
    file = os.path.join(anal_dir, 'eWES', subDir, sample, 'SNV', 'somatic', 'lofreq', sample + '.lofreq.cleaned.vcf')
    if not os.path.isfile(file) :
        print("lofreq: Intermediate file does not exist.")
    else :
        data = search_vcf(file, locus[0], locus[1], window)
        if data.shape[0] == 0 :
            print('lofreq: No matching data found.')
        else :
            print('[lofreq]')
            print(data.to_string(index=False))

    # freebayes
    file = os.path.join(anal_dir, 'eWES', subDir, sample, 'SNV', 'somatic', 'freebayes', sample + '.freebayes.cleaned.vcf')
    if not os.path.isfile(file) :
        print("freebayes: Intermediate file does not exist.")
    else :
        data = search_vcf(file, locus[0], locus[1], window)
        if data.shape[0] == 0 :
            print('freebayes: No matching data found.')
        else :
            print('[freebayes]')
            print(data.to_string(index=False))


