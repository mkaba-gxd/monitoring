import os
import sys
import datetime
import pandas as pd
import numpy as np
from pathlib import Path
from .func import *

use_column_info = ['SAMPLE_ID','seqDir','PATIENT_NO','GENDER','BIRTH_DATE','AGE','SAMPLING_DATE','Clinician','OCCURRED_ORGAN','DIAGNOSIS_NAME','PRJ_TYPE','Cohort','Timepoint','Institution','BIOPSY_OR_SURGERY']
use_column_data = ['seqDir','SAMPLE_ID','CHROM','POS','REF','ALT','AF','DP','ALT_AD','SYMBOL','HGVSc','HGVSp','Clinvar_CLNSIG','ONCOKB_ONCOGENICITY','FILTER_AF_gnomADg','FILTER_AF_gnomADe','FILTER_AF_tommo','FILTER_NONCODING','FILTER_SPLICING','FILTER_SYNONYMOUS','FILTER_BLACKLIST']
use_column_detail = ['seqDir','SAMPLE_ID','CHROM','POS','REF','ALT','AF','DP','SYMBOL','HGVSc','HGVSp','Clinvar_CLNSIG','ONCOKB_ONCOGENICITY','REF_mutect2','ALT_mutect2','AF_mutect2','DP_mutect2','REF_lofreq','ALT_lofreq','AF_lofreq','DP_lofreq','REF_freebayes','ALT_freebayes','AF_freebayes','DP_freebayes','FILTER_ANNO','FILTER']
use_column_uniq = ['CHROM','POS','SYMBOL','FILTER_ANNO','FILTER']

def sampleid_query() :
    query = f"""
    SELECT tesh.run_id, tesh.fc_id, concat(tesh.equip_side, tesh.fc_id) AS sub_name, gp.SAMPLE_ID, gp.PATIENT_NO, gp.GENDER, gp.BIRTH_DATE, gp.AGE, gp.SAMPLING_DATE, gp.PI_NAME AS Clinician, gp.OCCURRED_ORGAN, gp.DIAGNOSIS_NAME, gp.PRJ_TYPE, tol.cohort AS Cohort, tol.timepoint AS Timepoint, cctm.CLINICAL_TRIAL_NAME AS Title, cpcm.PI_COMP_NAME AS Institution, tol.biopsy_or_surgery AS BIOPSY_OR_SURGERY
    FROM gxd.tb_expr_seq_header tesh
    INNER JOIN gxd.gc_qc_sample gqs
    ON tesh.run_id = gqs.run_id
    INNER JOIN gxd.gc_project gp
    ON gqs.SAMPLE_ID = gp.SAMPLE_ID
    INNER JOIN gxd.gc_history_log ghl
    ON gqs.SAMPLE_ID = ghl.SAMPLE_ID
    AND ghl.idx = (SELECT MAX(idx) FROM gc_history_log WHERE SAMPLE_ID = gqs.SAMPLE_ID)
    INNER JOIN gxd.tb_order_line tol
    ON tol.sample_id = gp.SAMPLE_ID
    LEFT OUTER JOIN gxd.tb_order_header toh
    ON tol.order_header_id = toh.order_header_id
    LEFT OUTER JOIN gxd.cm_pi_company_mst cpcm
    ON toh.pi_comp = cpcm.PI_COMP_ID
    LEFT OUTER JOIN cm_pi_company_clinical_trial_map cpcctm
    ON toh.pi_comp =  cpcctm.PI_COMP_ID
    LEFT OUTER JOIN gxd.cm_clinical_trial_mst cctm
    ON cpcctm.CLINICAL_TRIAL_ID = cctm.CLINICAL_TRIAL_ID
    WHERE gp.PRJ_TYPE = 'EWES' AND ghl.ANAL_STATUS ='102'
    """
    return query

def convert_af(x):

    if isinstance(x, str) and "," in x:
        return x

    try:
        return float(x)
    except:
        return x

def run_intermediate(args):

    flowcellid = args.flowcellid
    spread = args.spread
    frequency = args.frequency
    directory = args.directory
    outdir = args.outdir

    if not os.path.exists(outdir): os.makedirs(outdir)
    now = datetime.datetime.now()
    outfile = os.path.join(outdir, now.strftime("%Y%m%d%H%M") + '.3tools.xlsx')

    if os.path.isfile(outfile) :
        try :
            os.remove(outfile)
        except PermissionError:
            print('Failed to delete existing file (permission issue)')
            sys.exit()
        except OSError as e:
            print('OS error during file deletion')
            sys.exit()
        except Exception as e:
            sys.exit({e})

    df_info = getinfo(sampleid_query())
    if df_info.shape[0] == 0 : init()

    if flowcellid is not None :
        df_info = df_info[ df_info['fc_id']==flowcellid ]
        if df_info.shape[0] == 0 : init()

    df_info = df_info[ df_info['Title']=='MONSTAR-SCREEN-3' ]
    if df_info.shape[0] == 0 : init()

    df_info = df_info[~df_info['SAMPLE_ID'].str.contains('_PCE_|_NCE_|_PCT_|_NCT_', regex=True, na=False)]
    if df_info.shape[0] == 0 : init()

    df_info['PRJ_TYPE'] = df_info['PRJ_TYPE'].str.replace('EWES',"eWES")
    df_info['Institution'] = df_info['Institution'].str.replace('　','')

    uniq_info = fcDir_table(df_info, directory)
    if uniq_info.shape[0] == 0: init()

    df_info = pd.merge(df_info, uniq_info, on=['sub_name','PRJ_TYPE'])
    df_info = df_info[use_column_info].drop_duplicates()
    df_info = df_info.sort_values(['SAMPLE_ID']).reset_index(drop=True)

    df_spread = None
    df_result = None
    df_detail = pd.DataFrame(columns=use_column_detail)

    for i, item in df_info.iterrows() :

        tempfile = os.path.join(directory, 'eWES', item['seqDir'], item['SAMPLE_ID'], 'SNV', 'somatic', item['SAMPLE_ID'] + '.target.snv.marked.tsv')
        if not os.path.isfile(tempfile) : continue

        data = pd.read_csv(tempfile, sep="\t", low_memory=False, dtype=str)
        data['seqDir'] = item['seqDir']
        data['SAMPLE_ID'] = item['SAMPLE_ID']
        data = data[use_column_data].drop_duplicates()
        data["HGVSc"] = data["HGVSc"].str.split(":", expand=True)[1]
        data["HGVSp"] = data["HGVSp"].str.split(":", expand=True)[1]

        data["FILTER_ANNO"] = np.where(
        data["Clinvar_CLNSIG"].str.contains("Pathogenic|Likely_pathogenic", case=True, na=False) |
        data["ONCOKB_ONCOGENICITY"].str.contains("oncogenic|resistance", case=False, na=False), "PASS", "FAIL")

        data["FILTER"] = "FAIL"
        data.loc[
        (data["FILTER_AF_gnomADg"] == "PASS") &
        (data["FILTER_AF_gnomADe"] == "PASS") &
        (data["FILTER_AF_tommo"] == "PASS") &
        (data["FILTER_BLACKLIST"] == "PASS") &
        (data["FILTER_NONCODING"] == "PASS") &
        ((data["FILTER_SYNONYMOUS"] == "PASS") | ((data["FILTER_SYNONYMOUS"] == "FAIL") & (data["FILTER_SPLICING"] == "FAIL"))) &
        True, "FILTER"] = "PASS"

        data = data.sort_values('SYMBOL').reset_index(drop=True)
        data['POS'] = pd.to_numeric(data['POS'], errors='coerce')
        data['AF'] = pd.to_numeric(data['AF'], errors='coerce')
        data['DP'] = pd.to_numeric(data['DP'], errors='coerce')
        data['ALT_AD'] = pd.to_numeric(data['ALT_AD'], errors='coerce')

        if spread or frequency :
            if df_spread is None :
                df_spread = data
            else :
                df_spread = pd.concat([df_spread, data], axis=0)

        df_filt = data[ data['AF'] >= 0.05 ]
        df_filt = df_filt[ df_filt['DP'] >= 20 ]
        df_filt = df_filt[ df_filt['ALT_AD'] >= 2 ]
        df_filt = df_filt[ df_filt['FILTER']=='PASS' ]
        if df_filt.shape[0] == 0 : continue

        if df_result is None :
            df_result = df_filt
        else :
            df_result = pd.concat([df_result, df_filt], axis=0)

        df_report = df_filt[ df_filt['FILTER_ANNO']=='PASS' ].reset_index(drop=True)
        if df_report.shape[0] == 0 : continue

        # mutect2
        file = os.path.join(directory, 'eWES', item['seqDir'], item['SAMPLE_ID'], 'SNV', 'somatic', 'mutect2', item['SAMPLE_ID'] + '.mutect2.cleaned.vcf')
        if not os.path.isfile(file) :
            dt_mt = None
        else :
            dt_mt = load_vcf(file, True)

        # lofreq
        file = os.path.join(directory, 'eWES', item['seqDir'], item['SAMPLE_ID'], 'SNV', 'somatic', 'lofreq', item['SAMPLE_ID'] + '.lofreq.cleaned.vcf')
        if not os.path.isfile(file) :
            dt_lr = None
        else :
            dt_lr = load_vcf(file, False)

        # freebayes
        file = os.path.join(directory, 'eWES', item['seqDir'], item['SAMPLE_ID'], 'SNV', 'somatic', 'freebayes', item['SAMPLE_ID'] + '.freebayes.cleaned.vcf')
        if not os.path.isfile(file) :
            dt_fb = None
        else :
            dt_fb = load_vcf(file, True)

        for j, value in df_report.iterrows() :
            val_mt = search_vcf(dt_mt, value['CHROM'], value['POS'], value['REF'], value['ALT'])
            val_lr = search_vcf(dt_lr, value['CHROM'], value['POS'], value['REF'], value['ALT'])
            val_fb = search_vcf(dt_fb, value['CHROM'], value['POS'], value['REF'], value['ALT'])

            new_row = [item['seqDir'], item['SAMPLE_ID']]
            new_row.extend(value[['CHROM','POS','REF','ALT','AF','DP','SYMBOL','HGVSc','HGVSp','Clinvar_CLNSIG','ONCOKB_ONCOGENICITY']].tolist())
            new_row.extend(val_mt)
            new_row.extend(val_lr)
            new_row.extend(val_fb)
            new_row.extend([value['FILTER_ANNO'],value['FILTER']])
            df_detail.loc[len(df_detail)] = new_row

    df_detail['AF_mutect2'] = df_detail['AF_mutect2'].apply(convert_af)
    df_detail['DP_mutect2'] = pd.to_numeric(df_detail['DP_mutect2'], errors='coerce')
    df_detail['AF_lofreq'] = pd.to_numeric(df_detail['AF_lofreq'], errors='coerce')
    df_detail['DP_lofreq'] = pd.to_numeric(df_detail['DP_lofreq'], errors='coerce')
    df_detail['DP_freebayes'] = pd.to_numeric(df_detail['DP_freebayes'], errors='coerce')
    df_detail = df_detail.sort_values(by=['seqDir','SAMPLE_ID','SYMBOL','HGVSc']).reset_index(drop=True)

    if frequency :
        df_freq = df_spread[ df_spread['FILTER']=='PASS' ].groupby(["CHROM", "POS", "SYMBOL"]).size().reset_index(name="count")
        df_freq = df_freq.sort_values(by="count", ascending=False)

    with pd.ExcelWriter(outfile, engine='openpyxl') as writer:
        df_info.to_excel(writer, sheet_name="index", index=False)
        if spread : df_spread.to_excel(writer, sheet_name="spread", index=False)
        if frequency : df_freq.to_excel(writer, sheet_name="frequency", index=False)
        df_detail.to_excel(writer, sheet_name="reported", index=False)


