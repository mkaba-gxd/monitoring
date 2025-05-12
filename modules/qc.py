import os
import sys
import datetime
import pandas as pd
from .func import *

class sql_code():
    def extracting_qc_eWES():
        query = f'''
        select DISTINCT l.sample_id, l.arrival_dt, COALESCE(l.ref_prep_id, l.prep_id) AS prep_id, l.customer_sample_id,
        -- sample prep
        ROUND(COALESCE(p1.qubit_conc,p2.qubit_conc),2)   AS DNA_QUBIT_CONC,
        ROUND(COALESCE(p1.final_amt,p2.final_amt),2)     AS FINAL_AMT,
        ROUND(COALESCE(p1.median_size,p2.median_size),2) AS MEDIAN_SIZE,
        ROUND(COALESCE(p1.DIN,p2.DIN),2)                 AS DIN,
        COALESCE(p1.grade,p2.grade)                      AS GRADE,
        -- library prep
        pre.fragment_size AS PRE_FRAGMENT_SIZE,
        pre.final_amt AS PRE_FINAL_AMT,
        pre.avg_size AS PRE_AVG_SIZE,
        pre.grade AS PRE_GRADE,
        post.mol AS MOL,
        post.avg_size AS POST_AVG_SIZE,
        post.grade AS POST_GRADE,
        q.FC_ID AS FC_ID,
        -- raw fastq
        b.FASTQ_TOTAL_READ_R1, -- raw fastq total read (before)
        b.FASTQ_GC_CONTENTS_R1, -- raw fastq GC contents (before)
        b.FASTQ_Q30_R1,        -- raw fastq Q30 (before)
        b.FASTQ_TOTAL_READ_R2, -- raw fastq total read (after)
        b.FASTQ_GC_CONTENTS_R2, -- raw fastq GC contents (after)
        b.FASTQ_Q30_R2,         -- raw fastq Q30 (after)
        -- bam
        b.BAM_MEAN_DEPTH_TARGET, -- bam mean depth target region
        b.BAM_MEAN_DEPTH AS BAM_MEAN_DEPTH_EXOME, -- bam mean depth exome region
        b.DEPTH_OF_COVERAGE, -- bam depth of coverage
        b.BAM_COVERAGE_100X, -- bam coverage 100x Rate
        b.BAM_ON_TARGET_RATE, -- bam on target rate
        b.BAM_DUP_RATE, -- bam duplication rate
        b.BAM_UNIFORM -- bam uniformity
        from gxd.tb_order_line l
        left join tb_expr_prep p1
        on l.ref_prep_id = p1.prep_id
        left join tb_expr_prep p2
        on l.prep_id = p2.prep_id
        left join tb_expr_pre_pcr pre
        on l.sample_id = pre.sample_id
        left join tb_expr_post_pcr post
        on l.sample_id = post.sample_id
        left join gc_qc_sample q
        on l.sample_id = q.sample_id
        left join gc_qc_bi b
        on l.sample_id = b.sample_id
        and b.idx = (select max(idx) from gc_qc_bi where sample_id = l.sample_id)
        WHERE l.sample_id regexp '^CD_'
        ORDER BY l.sample_id
        '''
        return query

    def extracting_qc_WTS():
        query = f'''
        select DISTINCT l.sample_id, l.arrival_dt, COALESCE(l.ref_prep_id, l.prep_id) AS prep_id, l.customer_sample_id,
        ROUND(COALESCE(p1.uv_conc,p2.uv_conc),2)       AS RNA_UV_CONC,
        ROUND(COALESCE(p1.final_amt,p2.final_amt),2)   AS FINAL_AMT,
        ROUND(COALESCE(p1.dv200,p2.dv200),2)           AS DV200,
        ROUND(COALESCE(p1.RIN,p2.RIN),2)               AS RIN,
        ROUND(COALESCE(p1.ratio_260_280,p2.ratio_260_280),2) AS RATIO260_280,
        COALESCE(p1.grade,p2.grade)                    AS GRADE,
        IFNULL(pre.final_amt,(ROUND((pre.qubit_conc * pre.final_vol)/1000,3))) AS PRE_FINAL_AMT,
        pre.avg_size AS PRE_AVG_SIZE,
        pre.grade AS PRE_GRADE,
        post.mol AS POST_MOL,
        post.avg_size AS POST_AVG_SIZE,
        post.grade AS POST_GRADE,
        q.FC_ID AS FC_ID,
        -- raw fastq
        b.FASTQ_TOTAL_READ_R1, -- raw fastq total read (before)
        b.FASTQ_GC_CONTENTS_R1, -- raw fastq GC contents (before)
        b.FASTQ_Q30_R1,        -- raw fastq Q30 (before)
        b.FASTQ_TOTAL_READ_R2, -- raw fastq total read (after)
        b.FASTQ_GC_CONTENTS_R2, -- raw fastq GC contents (after)
        b.FASTQ_Q30_R2,         -- raw fastq Q30 (after)
        -- bam
        b.BAM_UNIQUE_MAPPED_READS_PERCENT, -- bam unique mapped reads
        b.BAM_MULTIMAPPED_READS_PERCENT, -- bam multimapped reads
        b.BAM_UNMAPPED_READS_PERCENT, -- bam unmapped reads
        b.BAM_CHIMERIC_READS_PERCENT -- bam chimeric reads
        from gxd.tb_order_line l
        left join tb_expr_prep p1
        on l.ref_prep_id = p1.prep_id
        left join tb_expr_prep p2
        on l.prep_id = p2.prep_id
        left join tb_expr_pre_pcr pre
        on l.sample_id = pre.sample_id
        left join tb_expr_post_pcr post
        on l.sample_id = post.sample_id
        left join gc_qc_sample q
        on l.sample_id = q.sample_id
        left join gc_qc_bi b
        on l.sample_id = b.sample_id
        and b.idx = (select max(idx) from gc_qc_bi where sample_id = l.sample_id)
        WHERE l.sample_id regexp '^CR_'
        ORDER BY l.sample_id
        '''
        return query

def run_qc(args):

    out_file = args.output
    if out_file == "":
        now = datetime.datetime.now()
        out_file = '/data1/work/monitoring/QC/' + now.strftime("%Y%m%d%H%M") + '.xlsx'

    out_file = os.path.abspath(out_file)
    dir_path = os.path.dirname(out_file)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    if not out_file.lower().endswith(".xlsx"):
        out_file += ".xlsx"

    for testtype in ['eWES','WTS'] :
        data = getinfo(getattr(sql_code, f"extracting_qc_{testtype}")())
        try:
            with pd.ExcelWriter(out_file, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer :
                data.to_excel(writer, sheet_name=testtype, index=False)
        except FileNotFoundError:
            with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
                data.to_excel(writer, sheet_name=testtype, index=False)


