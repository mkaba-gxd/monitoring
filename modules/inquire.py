import os
import sys
import re
import glob
import warnings
import datetime
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
from .func import *

def splitLocus(locus) :

    pattern = re.compile(r"[A-Za-z0-9]+:[0-9]+-[0-9]+")
    valid_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}

    if not pattern.fullmatch(locus) :
        init("wrong locus")

    chrom = locus.split(':')[0]
    if chrom not in valid_chroms :
        init("wrong locus")

    start = locus.split(':')[1].split('-')[0]
    if start.isdigit() :
        start = int(start)
    else:
        init("wrong locus")

    stop  = locus.split(':')[1].split('-')[1]
    if stop.isdigit() :
        stop = int(stop)
    else :
        init("wrong locus")

    if stop <= start :
        init("wrong locus")


def create_inq(bam, forward: Path, tempDir, files=None, locus=None, sort=False) :

    os.makedirs(tempDir, exist_ok=True)
    bam_name = os.path.splitext(os.path.basename(bam))[0]
    tmp_name = os.path.basename(tempDir)

    cmd = ''
    idx_flag = True
    mv_flag  = True

    if locus is None :
        if sort :
            if bam.endswith('.bam') :
                cmd = f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools sort -@ 12 -o {tempDir}/{bam_name}.bam {bam} && "
            elif bam.endswith('.sam') :
                cmd = f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools view -@ 12 -bh {bam} > {tempDir}/{bam_name}.tmp.bam && "
                cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools sort -@ 12 -o {tempDir}/{bam_name}.bam {tempDir}/{bam_name}.tmp.bam && "
                cmd += f"rm -rf {tempDir}/{bam_name}.tmp.bam && "
        elif os.path.isfile(bam + '.bai') :
            os.symlink(bam, os.path.join(tempDir, os.path.basename(bam)))
            os.symlink(bam + '.bai', os.path.join(tempDir, os.path.basename(bam) + '.bai'))
            idx_flag = False
            mv_flag  = False
        else :
            os.symlink(bam, os.path.join(tempDir, bam_name + '.bam'))
            mv_flag  = False

    else :

        splitLocus(locus)

        if sort :
            if bam.endswith('.bam') :
                cmd = f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools sort -@ 12 -o {tempDir}/{bam_name}.tmp.bam {bam} && "
                cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools index -@ 12 {tempDir}/{bam_name}.tmp.bam && "
                cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools view -@ 12 -bh {tempDir}/{bam_name}.tmp.bam {locus} > {tempDir}/{bam_name}.bam && "
                cmd += f"rm -rf {tempDir}/{bam_name}.tmp.bam {tempDir}/{bam_name}.tmp.bam.bai && "
            elif bam.endswith('.sam') :
                cmd = f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools view -@ 12 -bh {bam} > {tempDir}/{bam_name}.tmp.bam && "
                cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools sort -@ 12 -o {tempDir}/{bam_name}.sort.bam {tempDir}/{bam_name}.tmp.bam && "
                cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools index -@ 12 -bh {tempDir}/{bam_name}.sort.bam && "
                cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools view -@ 12 -bh {tempDir}/{bam_name}.sort.bam {locus} > {tempDir}/{bam_name}.bam && "
                cmd += f"rm -rf {tempDir}/{bam_name}.tmp.bam {tempDir}/{bam_name}.sort.bam {tempDir}/{bam_name}.sort.bam.idx && "

        elif os.path.isfile(bam + '.bai') :
            cmd = f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools view -@ 12 -bh {bam} {locus} > {tempDir}/{bam_name}.bam && "
        else :
            os.symlink(bam, os.path.join(tempDir, bam_name + '.tmp.bam'))
            cmd = f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools index -@ 12 {tempDir}/{bam_name}.tmp.bam &&"
            cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools view -@ 12 -bh {tempDir}/{bam_name}.tmp.bam {locus} > {tempDir}/{bam_name}.bam && "
            cmd += f"rm -rf {tempDir}/{bam_name}.tmp.bam {tempDir}/{bam_name}.tmp.bam.bai && "

    if idx_flag :
        cmd += f"singularity exec --bind /data1 /data1/GxD_eWES/Pipeline/containers/samtools.sif samtools index -@ 12 {tempDir}/{bam_name}.bam && "

    if not files is None :
        for sendFile in files :
            if os.path.isfile(sendFile) :
                cmd += f"rsync -azruL {sendFile} {tempDir}/ && "
            else :
                print("File not found: " + sendFile + "; skip")

    if mv_flag :
        cmd += f"mv {tempDir}/* {forward}/ && "
    else :
        cmd += f"rsync -azruL {tempDir}/* {forward}/ && "

    cmd += f"rm -rf {tempDir} "

    if os.path.isdir(forward) : 
        choice = prompt_choice("The output directory exists. Do you want to delete its contents? (yes[Y]/no[N]): ", ['yes', 'y', 'no', 'n'])
        if choice in ['yes', 'y']:
            shutil.rmtree(forward)

    os.makedirs(forward, exist_ok=True)

    qsubCmd = f"/data1/apps/sge/bin/lx-amd64/qsub -N INQ_{tmp_name} -q all.q -pe smp 12 -o /dev/null -e /dev/null << EOF\n{cmd}\nEOF"
#    print(qsubCmd)
    os.system(qsubCmd)

def run_inquire(args) :

    sample = args.sample
    item = args.item
    locus = args.locus
    outdir = args.outdir
    anal_dir = args.directory

    tempDir = Path(os.path.abspath(__file__)).parent.parent / "tmp"
    now = datetime.datetime.now()
    outdir = os.path.join(outdir, now.strftime("%Y%m%d"), sample)
    tempDir = os.path.join(tempDir, now.strftime("%d%H%M%S") + str(now.microsecond))

    if item == "snv" :
        subDir = batch(sample, anal_dir, 'eWES')
    elif item == "fusion" or item == "splice" :
        subDir = batch(sample, anal_dir, 'WTS')
    else :
        init("Wrong argument.")

    if subDir == "Incorrect" :
        init("No registration in database")

    if item == "snv" :
        subDir = os.path.join(anal_dir, 'eWES', subDir)
    else :
        subDir = os.path.join(anal_dir, 'WTS', subDir)

    if item == "snv" :
        bam = os.path.join(subDir, sample, 'Preprocessing','align',sample + '.tumour.recaled.bam')
        if not os.path.isfile(bam) :
            init('bam file does not exist: ' + bam)

        create_inq(bam=bam, locus=locus, sort=False, forward=outdir, tempDir=tempDir)

    elif item == "fusion" :
        bam = os.path.join(subDir, sample, 'Fusion','STAR-Fusion','STAR_align_starfu', sample + '.star-fusion.Aligned.out.bam')
        if not os.path.isfile(bam) :
            bam = os.path.join(subDir, sample, 'Fusion','STAR-Fusion','STAR_align_starfu', sample + '.star-fusion.Aligned.out.sam')
        if not os.path.isfile(bam) :
            init('sam/bam file does not exist: ' + bam)

        files = [ os.path.join(subDir, sample, 'Fusion', sample + '.fusion.filtered.tsv'), os.path.join(subDir, sample, 'Fusion', 'Arriba', sample + '.fusions.tsv'), os.path.join(subDir, sample, 'Fusion', 'STAR-Fusion', 'star-fusion.fusion_predictions.abridged.coding_effect.tsv')]

        create_inq(bam=bam, files=files, locus=locus, sort=True, forward=outdir, tempDir=tempDir)
 
    elif item == "splice" :
        bam = os.path.join(subDir, sample, 'Expression','STAR_align_exp',sample + '.Aligned.sortedByCoord.out.bam')
        if not os.path.isfile(bam) :
            init('bam file does not exist: ' + bam)

        create_inq(bam=bam, locus=locus, sort=False, forward=outdir, tempDir=tempDir)


