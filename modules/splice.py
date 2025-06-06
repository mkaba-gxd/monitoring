import os
import sys
import re
import glob
import math
import itertools
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from .func import *

#bed_f = os.path.join(Path(os.path.abspath(__file__)).parent, 'AR_MET.bed')
#target = os.path.join(Path(os.path.abspath(__file__)).parent, 'AR_MET_target.txt')
rscript = os.path.join(Path(os.path.abspath(__file__)).parent, 'covFigure.R')

#if not os.path.isfile(bed_f) : init(bed_f + ' file not found.')
#if not os.path.isfile(target) : init(target + ' file not found.')
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
    dpCmd = rsCmd = ''
    mvCmd = f"cd {outdir} && "
    for c, p in BED.items():
        dpCmd = dpCmd + f"samtools depth -@ 12 -a -Q 255 -b {p} {bam} > {sample}.{c}.depth.txt && "
    for c, p in TARGET.items():
        rsCmd = rsCmd + f"Rscript --no-save --slave {rscript} {sample}.{c}.depth.txt {p} {sample} && "
    cmd = mvCmd + dpCmd + rsCmd.rstrip('&& ')
    os.system(cmd)
    
