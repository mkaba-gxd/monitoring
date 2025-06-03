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

'''
def batch(sample, anal_dir):

    tbl = getinfo(subname_query(sample))
    if tbl.shape[0] == 0 : init("Unregistered sample ID.")
    subname = tbl.sub_name[0]
    anal_type = tbl.PRJ_TYPE[0]
    if anal_type != 'WTS' : init("Analysis type is not WTS")
    anal_dir = Path(os.path.join(anal_dir, anal_type))
    fcDirs = [fcDir for fcDir in anal_dir.iterdir() if fcDir.name.endswith(subname)]
    fcDirs.sort()

    return os.path.basename(fcDirs[-1])

def subname_query(sample):
    query = f"""
    SELECT concat(tesh.equip_side, tesh.fc_id) AS sub_name, gp.PRJ_TYPE
    FROM gxd.tb_expr_seq_header tesh
    INNER JOIN gxd.gc_qc_sample gqs
    ON tesh.run_id = gqs.run_id
    INNER JOIN gxd.gc_project gp
    ON gqs.SAMPLE_ID = gp.SAMPLE_ID
    INNER JOIN gxd.gc_history_log ghl
    ON gqs.SAMPLE_ID = ghl.SAMPLE_ID
    AND ghl.idx = (SELECT MAX(idx) FROM gc_history_log WHERE SAMPLE_ID = gqs.SAMPLE_ID)
    WHERE ghl.SAMPLE_ID = '{sample}'
    """
    return query
'''

def read_line(file_path):
    if not os.path.isfile(file_path) : return 0
    cmd = f"wc -l {file_path}"
    c = int(subprocess.check_output(cmd.split()).decode().split()[0]) / 2
    return c

def run_seqr(args):

    sample = args.sample
    verbose = args.verbose
    anal_dir = args.analysis_dir

    subDir = batch(sample, anal_dir)
    if subDir is None : init("No registration in database")

    chim_dir = os.path.join(anal_dir, 'WTS', subDir, sample, 'Fusion', 'STAR-SEQR', f"{sample}_STAR-SEQR", 'chim_transcripts')
    if not os.path.isdir(chim_dir) : init("chim_transcripts folder not yet created: " + chim_dir)

    files_file = glob.glob(chim_dir + '/transcripts-fusion*.fa')
    clean_file = [ os.path.split(x)[1] for x in files_file ]
    clean_file = [ re.sub('transcripts-fusion-','',x) for x in clean_file ]
#    clean_file = [ re.sub('transcripts-left-','',x) for x in clean_file ]
#    clean_file = [ re.sub('transcripts-right-','',x) for x in clean_file ]

    brk1 = pd.DataFrame({'brk':[ '_'.join(x.split('_')[0:3]) for x in clean_file ], 'file_path':files_file})
    brk2 = pd.DataFrame({'brk':[ '_'.join(x.split('_')[3:6]) for x in clean_file ], 'file_path':files_file})
    all_brk_map = pd.concat([brk1,brk2]).reset_index().drop('index', axis=1)
    shared_brks = all_brk_map[all_brk_map.duplicated(subset='brk')]['brk'].unique()

    total = 0
    for brk in shared_brks :
        part_brk_map = all_brk_map[all_brk_map['brk']==brk]
        lines = [ read_line(f) for f in part_brk_map['file_path'] ]
        comb = list(itertools.combinations(lines, 2))
        total += int(sum([ math.prod(x) for x in comb ]))
        if verbose :
            temp_brk_map = pd.DataFrame({'index':[ re.sub(".fa$",'',re.sub('transcripts-fusion-','',os.path.basename(f))) for f in part_brk_map['file_path'] ], 'counts':lines})
            print(temp_brk_map)

    print('Combinations of all sequences: ' + f'{total:,}')
#    print(sample + "\t" + f'{total:,}')

    
