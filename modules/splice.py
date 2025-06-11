import os
import sys
import re
import glob
import math
import warnings
import itertools
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from .func import *

rscript = os.path.join(Path(os.path.abspath(__file__)).parent, 'covFigure.R')
if not os.path.isfile(rscript) : init(rscript + ' file not found.')

LIST = ['EGFR','MET','AR']

def set_file_path(category, filetype):

    FILE_PATH = {}
    for c in category :
        path = os.path.join(Path(os.path.abspath(__file__)).parent.parent, 'reference', c + '.' + filetype)
        if os.path.isfile(path) : 
            FILE_PATH[c] = path

    return FILE_PATH


def run_splice(args) :

    sample = args.sample
    category = args.category
    outdir = args.outdir
    anal_dir = args.analysis_dir

    subDir = batch(sample, anal_dir)
    if subDir is None : init("No registration in database")

    inc = [item for item in category if item in LIST ]
    exc = [item for item in category if item not in LIST ]

    if len(exc) > 0 :
        print(','.join(exc) + ' is not covered.')

    if len(inc) == 0 :
        init('No category. Exit.')

    bam = os.path.join(anal_dir, 'WTS', subDir, sample, 'Expression', 'STAR_align_exp', sample + '.Aligned.sortedByCoord.out.bam')
    if not os.path.isfile(bam) : init("bam file not yet created: " + bam)

    BED = set_file_path(category,'bed')
    TARGET = set_file_path(category,'txt')
    if len(BED) != len(category) or len(TARGET) != len(category) :
        init('Reference file not found. Abort.')

    os.makedirs(outdir, exist_ok=True)
    dpCmd = rsCmd = rmCmd = ''
    mvCmd = f"cd {outdir} && "
    for c, p in BED.items():
        dpCmd = dpCmd + f"samtools depth -@ 12 -a -Q 255 -b {p} {bam} > {sample}.{c}.depth.txt && "
    for c, p in TARGET.items():
        rsCmd = rsCmd + f"Rscript --no-save --slave {rscript} {sample}.{c}.depth.txt {p} {sample} && "
        rmCmd = rmCmd + f"rm {sample}.{c}.depth.txt && " 
    cmd = mvCmd + dpCmd + rsCmd + rmCmd.rstrip('&& ')
    os.system(cmd)

    result_file = os.path.join(anal_dir, 'WTS', subDir, sample, 'Alternative_splicing','ESDetector', sample + '.exon_skipped.annotated.tsv')
    use_column = ['spliceName','discordant_mates','canonical_reads','ratio','tpm_total','tpm_variant']
    out_df = pd.DataFrame(columns=use_column)

    try :
        anno = pd.read_csv(result_file, sep="\t")
        anno = anno[use_column]
        for c in category :
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                out_df = pd.concat([out_df, anno[anno['spliceName'].str.contains(c)]], axis=0, ignore_index=True)
    except Exception as e:
        init(e)
 
    if out_df.shape[0] == 0:
        init('No ESDetector results.')

    print(out_df)

