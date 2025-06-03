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

bed_f = os.path.join(Path(os.path.abspath(__file__)).parent, 'AR_MET.bed')
target = os.path.join(Path(os.path.abspath(__file__)).parent, 'AR_MET_target.txt')
rscript = os.path.join(Path(os.path.abspath(__file__)).parent, 'covFigure.R')

if not os.path.isfile(bed_f) : init(bed_f + ' file not found.')
if not os.path.isfile(target) : init(target + ' file not found.')
if not os.path.isfile(rscript) : init(rscript + ' file not found.')

def run_splice(args) :

    sample = args.sample
    outdir = args.outdir
    anal_dir = args.analysis_dir

    subDir = batch(sample, anal_dir)
    if subDir is None : init("No registration in database")

    bam = os.path.join(anal_dir, 'WTS', subDir, sample, 'Expression', 'STAR_align_exp', sample + '.Aligned.sortedByCoord.out.bam')
    if not os.path.isfile(bam) : init("bam file not yet created: " + bam)

    os.makedirs(outdir, exist_ok=True)
    mvCmd = f"cd {outdir} && "
    dpCmd = f"samtools depth -@ 12 -a -Q 255 -b {bed_f} {bam} > {sample}.depth.txt && "
    rsCmd = f"Rscript --no-save --slave {rscript} {sample}.depth.txt {target} {sample}"
    cmd = mvCmd + dpCmd + rsCmd
    os.system(cmd)


    
