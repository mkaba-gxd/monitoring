import pandas as pd
import datetime
import warnings
from .func import *

def remove_dup_list(lst):
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]

def m3_query() :
    query = f"""
    SELECT tesh.run_id, concat(tesh.equip_side, tesh.fc_id) AS sub_name, gp.PRJ_TYPE, gp.DIAGNOSIS_NAME, gp.PATIENT_NO, gp.SAMPLE_ID
    FROM gxd.tb_expr_seq_header tesh
    INNER JOIN gxd.gc_qc_sample gqs
    ON tesh.run_id = gqs.run_id
    INNER JOIN gxd.gc_project gp
    ON gqs.SAMPLE_ID = gp.SAMPLE_ID
    INNER JOIN gxd.gc_history_log ghl
    ON gqs.SAMPLE_ID = ghl.SAMPLE_ID
    AND ghl.idx = (SELECT MAX(idx) FROM gc_history_log WHERE SAMPLE_ID = gqs.SAMPLE_ID)
    WHERE gp.PRJ_TYPE = 'EWES' AND ghl.ANAL_STATUS ='102'
    """
    return query

def run_cnv(args):

    genes = [x.strip() for x in args.genes.split(',') if not x.strip() == '']
    directory = args.directory
    outdir = args.outdir
    exclusion = [x.strip() for x in args.exclusion.split(',') if not x.strip() == '']
    exclusion = rmdup_list(exclusion)

    now = datetime.datetime.now()
    out_file = os.path.join(outdir, now.strftime("%Y%m%d%H%M") + '.xlsx')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    genes = remove_dup_list(genes)
    if len(genes) == 0: init('Gene name input value error')

    df_info = getinfo(m3_query())
    if df_info.shape[0] == 0 : init()
    df_info = df_info[ df_info['PATIENT_NO'].str.match( r'^M3\d{7}$', na=False) ]
    if df_info.shape[0] == 0 : init()

    if len(exclusion) > 0:
        print ("exclusion sample:" + ",".join(exclusion))
        df_info = df_info[ ~df_info['SAMPLE_ID'].isin(exclusion)]
        if df_info.shape[0] == 0 : init("No corresponding sample IDs.")

    df_info['PRJ_TYPE'] = df_info['PRJ_TYPE'].str.replace('EWES',"eWES")
    df_info = df_info[ df_info['PRJ_TYPE'] == "eWES" ]
    if df_info.shape[0] == 0 : init("Test type error: no sample ID corresponds.")

    df_info = df_info[~df_info['SAMPLE_ID'].str.contains('_PCE_|_NCE_|_PCT_|_NCT_', regex=True, na=False)]
    if df_info.shape[0] == 0 : init("No clinical specimens match the criteria.")

    uniq_info = fcDir_table(df_info, directory)
    if uniq_info.shape[0] == 0: init()

    df_info = pd.merge(df_info, uniq_info, on=['sub_name','PRJ_TYPE'])

    merge_data = None

    for i, item in df_info.iterrows() :
        anal_dir = os.path.join(directory,"eWES",item['seqDir'],item['SAMPLE_ID'],'CNV')
        cnv_file = os.path.join(anal_dir, item['SAMPLE_ID'] + '.cnv.marked.tsv')

        try :
            df = pd.read_csv(cnv_file, sep="\t", header=0, low_memory=False)
            df = df[['Gene_name','CHROM','START','END','gene.mean.CN']]
            
            filt = df[ df['Gene_name'].isin(genes) ].drop_duplicates().copy()
            miss = list(set(genes) - set(filt['Gene_name']))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=FutureWarning)

                if miss :
                    na_rows = pd.DataFrame({'Gene_name': miss})
                    for col in df.columns[ df.columns != 'Gene_name' ]:
                        na_rows[col] = pd.NA
                    df = pd.concat([filt, na_rows], ignore_index=True)
                else :
                    df = filt

            df.insert(0, 'sample_id', item['SAMPLE_ID'])

            if merge_data is None :
                merge_data = df
            else:
                merge_data = pd.concat([merge_data, df], axis=0)

        except Exception as e:
            continue

    miss = []
    for g in genes :
        df_g = merge_data[ merge_data['Gene_name']==g ]
        if df_g.dropna(subset=['gene.mean.CN']).shape[0] == 0 :
            miss.append(g)
            continue

        df_g = df_g.sort_values(['sample_id']).reset_index(drop=True)
        try :
            with pd.ExcelWriter(out_file, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer :
                df_g.to_excel(writer, sheet_name=g, index=False)
        except FileNotFoundError:
            with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
                df_g.to_excel(writer, sheet_name=g, index=False)

    if len(miss) > 0 :
        print('[' + ','.join(miss) + '] does not exist in the reference.')

