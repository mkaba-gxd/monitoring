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
use_columns = ["CHROM", "POS", "REF", "ALT","FORMAT","INFO"]
out_columns = ["CHROM", "POS", "REF", "ALT","AF","DP"]

def load_vcf(vcf_path, d_type):

    with open(vcf_path) as f:
        skip_rows = sum(1 for line in f if line.startswith('#'))

    if d_type :
        data = pd.read_csv(vcf_path, sep='\t', skiprows=skip_rows, header=None, index_col=False, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER','ATTR', 'FORMAT','INFO'])
    else :
        data = pd.read_csv(vcf_path, sep='\t', skiprows=skip_rows, header=None, index_col=False, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
        data['FORMAT'] = None

    data = data[use_columns].copy()

    return data

def search_vcf(vcf_path, d_type, chrom, position, window=0):

    data = load_vcf(vcf_path, d_type)
    data = data[ data['CHROM']==chrom ]
    if window == 0:
        data = data[data['POS']==int(position)].reset_index(drop=True)
    else :
        data = data[(int(position)-window<=data['POS']) & (data['POS']<=int(position)+window)].reset_index(drop=True)

    if data.shape[0] == 0 :
        return data

    return add_column(data)

def add_column(data) :

    data['AF'] = None
    data['DP'] = None

    for i, item in data.iterrows() :
        if item['FORMAT'] is None :
            info_dict = dict(x.split("=", 1) for x in item['INFO'].split(";") if "=" in x)
        else :
            info_keys = item['FORMAT'].split(":")
            info_value = item['INFO'].split(":")
            info_dict = {key: val for key, val in zip(info_keys, info_value)}

        data.loc[i, 'AF'] = info_dict.get("AF")
        data.loc[i, 'DP'] = info_dict.get("DP")

    return data[out_columns]

def search_vcf_mt2(vcf_path, chrom, position, window=0):

    with open(vcf_path) as f:
        skip_rows = sum(1 for line in f if line.startswith('#'))

    data = pd.read_csv(vcf_path, sep='\t', skiprows=skip_rows, header=None, index_col=False, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER','ATTR', 'FORMAT','INFO'])

    data = data[ data['CHROM']==chrom ]
    if window == 0:
        data = data[data['POS']==int(position)].reset_index(drop=True)
    else :
        data = data[(int(position)-window<=data['POS']) & (data['POS']<=int(position)+window)].reset_index(drop=True)

    if data.shape[0] == 0 :
        return data

    data['AF'] = None
    data['DP'] = None
    data['RPA'] = None

    for i, item in data.iterrows() :
        if item['FORMAT'] is None :
            info_dict = dict(x.split("=", 1) for x in item['INFO'].split(";") if "=" in x)
        else :
            info_keys = item['FORMAT'].split(":")
            info_value = item['INFO'].split(":")
            info_dict = {key: val for key, val in zip(info_keys, info_value)}

        data.loc[i, 'AF'] = info_dict.get("AF")
        data.loc[i, 'DP'] = info_dict.get("DP")

        match = re.search(r'RPA=([\d,]+)', item['ATTR'])
        if match:
            data.loc[i,'RPA'] = match.group(1)

    return data[["CHROM", "POS", "REF", "ALT","AF","DP","FILTER","RPA"]]


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
    if subDir == "Incorrect" : init("No registration in database")

    # mutect2
    file = os.path.join(anal_dir, 'eWES', subDir, sample, 'SNV', 'somatic', 'mutect2', sample + '.mutect2.cleaned.vcf')
    if not os.path.isfile(file) :
        print("mutect2: Intermediate file does not exist.")
        mt2_flag = True
    else :
        data = search_vcf(file, True, locus[0], locus[1], window)
        if data.shape[0] == 0 :
            print('mutect2: No matching data found.')
            mt2_flag = True
        else :
            print('[mutect2]')
            print(data.to_string(index=False))
            mt2_flag = False

    # lofreq
    file = os.path.join(anal_dir, 'eWES', subDir, sample, 'SNV', 'somatic', 'lofreq', sample + '.lofreq.cleaned.vcf')
    if not os.path.isfile(file) :
        print("lofreq: Intermediate file does not exist.")
    else :
        data = search_vcf(file, False, locus[0], locus[1], window)
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
        data = search_vcf(file, True, locus[0], locus[1], window)
        if data.shape[0] == 0 :
            print('freebayes: No matching data found.')
        else :
            print('[freebayes]')
            print(data.to_string(index=False))

    # If Mutect2 fails to detect
    if mt2_flag :
        file = os.path.join(anal_dir, 'eWES', subDir, sample, 'SNV', 'somatic', 'mutect2', sample + '.mutect2.filtered.vcf')
        if not os.path.isfile(file) : sys.exit(1)

        data = search_vcf_mt2(file, locus[0], locus[1], window)
        if data.shape[0] == 0 : sys.exit(1)
        print("\n[mutect2(Before filtering)]")
        print(data.to_string(index=False))


