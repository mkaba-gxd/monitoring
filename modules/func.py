import os
import sys
import pymysql
import warnings
import pandas as pd
from itertools import product
from pathlib import Path

def getinfo(comm):

    try :
        connection = pymysql.connect(host="192.168.9.100", user="gxd_pipeline", password="gw!2341234", database="gxd")
    except Exception as e:
        sys.exit({e})

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        db_tbl = pd.read_sql(comm, connection)

    return db_tbl

def SelectData(fc_id):
    query = f"""
    SELECT tesh.run_id, concat(tesh.equip_side, tesh.fc_id) AS sub_name, gp.PRJ_TYPE, gp.DIAGNOSIS_NAME, gp.PATH_NO, gp.SAMPLING_DATE, ghl.ANAL_STATUS, gp.SAMPLE_ID
    FROM gxd.tb_expr_seq_header tesh
    INNER JOIN gxd.gc_qc_sample gqs
    ON tesh.run_id = gqs.run_id
    INNER JOIN gxd.gc_project gp
    ON gqs.SAMPLE_ID = gp.SAMPLE_ID
    INNER JOIN gxd.gc_history_log ghl
    ON gqs.SAMPLE_ID = ghl.SAMPLE_ID
    AND ghl.idx = (SELECT MAX(idx) FROM gc_history_log WHERE SAMPLE_ID = gqs.SAMPLE_ID)
    WHERE tesh.fc_id = '{fc_id}'
    """
    return query

def SelectInfo(fc_id):
    query = f"""
    SELECT concat(tesh.equip_side, tesh.fc_id) AS sub_name, ghl.SAMPLE_ID, tol.timepoint, gp.PATIENT_NO, gp.PRJ_TYPE, gp.DIAGNOSIS_NAME, cctm.CLINICAL_TRIAL_NAME
    FROM gxd.tb_expr_seq_header tesh
    INNER JOIN gxd.gc_qc_sample gqs
    ON tesh.run_id = gqs.run_id
    INNER JOIN gxd.gc_project gp
    ON gqs.SAMPLE_ID = gp.SAMPLE_ID
    INNER JOIN gxd.tb_order_line tol
    ON tol.sample_ID = gp.SAMPLE_ID
    INNER JOIN gxd.gc_history_log ghl
    ON gqs.SAMPLE_ID = ghl.SAMPLE_ID
    AND ghl.idx = (SELECT MAX(idx) FROM gc_history_log WHERE SAMPLE_ID = gqs.SAMPLE_ID)
    LEFT OUTER JOIN gxd.tb_order_header toh
    ON tol.order_header_id = toh.order_header_id
    LEFT OUTER JOIN cm_pi_company_clinical_trial_map cpcctm
    ON toh.pi_comp =  cpcctm.PI_COMP_ID
    LEFT OUTER JOIN gxd.cm_clinical_trial_mst cctm
    ON cpcctm.CLINICAL_TRIAL_ID = cctm.CLINICAL_TRIAL_ID
    WHERE tesh.fc_id = '{fc_id}' AND cctm.CLINICAL_TRIAL_NAME = 'MONSTAR-SCREEN-3' AND ghl.ANAL_STATUS = '102'
    """
    return query

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

def fcDir_table(df, novaseqDir: Path):
    df = df[['sub_name','PRJ_TYPE']].drop_duplicates()
    df['seqDir'] = None

    for i, item in df.iterrows() :
        fcDir = Search_fcDir(item['sub_name'], Path(novaseqDir + '/' + item['PRJ_TYPE']))
        df.loc[i,'seqDir'] = fcDir

    df = df.dropna(subset=['seqDir'])

    return df

def Search_fcDir(batchID, novaseqDir : Path):

    fcDirs = [fcDir for fcDir in novaseqDir.iterdir() if fcDir.name.endswith(batchID)]
    fcDirs.sort()
    if len(fcDirs) != 1: return None

    return os.path.basename(fcDirs[-1])

def batch(sample, anal_dir, anal_type):

    tbl = getinfo(subname_query(sample))
    tbl['PRJ_TYPE'] = tbl['PRJ_TYPE'].str.replace('EWES',"eWES")
    if tbl.shape[0] == 0 : init("Unregistered sample ID.")
    if anal_type != tbl.PRJ_TYPE[0] : init("Analysis type is not " + anal_type)
    subname = tbl.sub_name[0]
    anal_dir = Path(os.path.join(anal_dir, anal_type))
    fcDirs = [fcDir for fcDir in anal_dir.iterdir() if fcDir.name.endswith(subname)]
    fcDirs.sort()

    return os.path.basename(fcDirs[-1])

def rmdup_list(lst):
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]

def expand_breakpoints(df):

    expanded_rows = []

    for _, row in df.iterrows():

        bp1_raw = str(row['breakpoint_1'])
        bp2_raw = str(row['breakpoint_2'])

        bp1_values = bp1_raw.split('|') if '|' in bp1_raw else [bp1_raw]
        bp2_values = bp2_raw.split('|') if '|' in bp2_raw else [bp2_raw]

        for bp1, bp2 in product(bp1_values, bp2_values):
            new_row = row.copy()
            new_row['breakpoint_1'] = bp1
            new_row['breakpoint_2'] = bp2
            expanded_rows.append(new_row)

    return pd.DataFrame(expanded_rows)

def load_vcf(vcf_path, d_type):

    with open(vcf_path) as f:
        skip_rows = sum(1 for line in f if line.startswith('#'))

    if d_type :
        data = pd.read_csv(vcf_path, sep='\t', skiprows=skip_rows, header=None, index_col=False, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER','ATTR', 'FORMAT','INFO'])
    else :
        data = pd.read_csv(vcf_path, sep='\t', skiprows=skip_rows, header=None, index_col=False, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
        data['FORMAT'] = None

    data = data[["CHROM", "POS", "REF", "ALT","FORMAT","INFO"]].copy()

    return data

def search_vcf(data, chrom, position, ref, alt):

    if data is None :
        return None, None, None, None

    data = data[ data['CHROM']==chrom ]
    data = data[data['POS']==int(position)]
    data = data[ data['REF']==ref]
    data = data[ data['ALT']==alt].reset_index(drop=True)

    if data.shape[0] == 0 :
        return None, None, None, None

    if data['FORMAT'][0] is None :
        info_dict = dict(x.split("=", 1) for x in data['INFO'][0].split(";") if "=" in x)
    else :
        info_keys = data['FORMAT'][0].split(":")
        info_value = data['INFO'][0].split(":")
        info_dict = {key: val for key, val in zip(info_keys, info_value)}

    return data['REF'][0], data['ALT'][0], info_dict.get("AF"), info_dict.get("DP")


def pic_value(file, column):
    try :
        df = pd.read_csv(file, sep="\t")
        if not column in df.columns :
            return '-'
        else :
            return df[column][0]
    except Exception as e:
        return '-'

def remove_files(FILES) :
    for file in FILES:
        if os.path.isfile(file):
            os.remove(file)

def init(msg="No matching data found.", parser=None):
    print(msg)
    if parser :
        parser.print_help()
    sys.exit(1)


