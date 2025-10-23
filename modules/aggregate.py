import os
import sys
import pandas as pd
import numpy as np
from itertools import product
from pathlib import Path
from .func import *

file_oncoKB = '/data1/GxD_eWES/reference/anno_db/oncoKB/v4.23/variant_table.tsv'

use_columns_ewes_1 = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','GT','REF_AD','ALT_AD','AF','DP','Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','MANE_SELECT','MANE_PLUS_CLINICAL','TSL','APPRIS','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','UNIPROT_ISOFORM','SOURCE','GENE_PHENO','SIFT','PolyPhen','DOMAINS','miRNA','HGVS_OFFSET','HGVSg','Clinvar','Clinvar_CLNSIG','Clinvar_ALLELEID','Clinvar_CLNDN','Clinvar_CLNHGVS','Clinvar_CLNVC','Clinvar_GENEINFO','Clinvar_MC','ONCOKB_ONCOGENICITY','ONCOKB_VARIANT_KEY','FILTER_DEPTH','FILTER_ALT_DEPTH','FILTER_VAF','FILTER_AF_gnomADg','FILTER_AF_gnomADe','FILTER_AF_tommo','FILTER_NONCODING','FILTER_SPLICING','FILTER_SYNONYMOUS','FILTER_BLACKLIST']
use_columns_ewes_2 = ['Gene_name','CHROM','START','END','TYPE','ONCOKB_VARIANT_CNV','gene.mean.CN','FILTER_CN','FILTER_ONCOKB','FILTER_GENES']

use_columns_wts_1 = ['samples','Out-of-Frame','OncoKB','cancer-related','gene1','gene2','chr1','breakpoint_1','chr2','breakpoint_2','max_split_cnt','max_span_cnt','sample_type','disease','tools','inferred_fusion_type','cancer_db_hits','fusion_IDs']
use_columns_wts_2 = ['gene1','gene2','strand1(gene/fusion)','strand2(gene/fusion)','breakpoint1','breakpoint2','discordant_mates','canonical_reads','ratio','tpm_total','tpm_variant']

def run_aggregate(args) :

    flowcellid = args.flowcellid
    directory = args.directory
    project_type = args.project_type
    outdir = args.outdir
    inclusion = [x.strip() for x in args.inclusion.split(',') if not x.strip() == '']
    exclusion = [x.strip() for x in args.exclusion.split(',') if not x.strip() == '']

    inclusion = rmdup_list(inclusion)
    exclusion = rmdup_list(exclusion)

    if len(inclusion) > 0 and len(exclusion) > 0 :
        init('ERROR: Inclusion and exclusion cannot be specified simultaneously.')

    var_tbl = pd.read_csv(file_oncoKB, sep="\t")
    var_tbl = var_tbl[var_tbl['VARIANT_GROUP_ID']=='SNV']
    var_tbl = var_tbl[['HUGO_SYMBOL','GRCH38_REFSEQ','GRCH38_ISOFORM']].drop_duplicates()
    var_tbl.columns = ['SYMBOL','ONCOKB_REFSEQ','ONCOKB_GENCODE']

    df_info = getinfo(SelectData(flowcellid))
    m3_info = getinfo(SelectInfo(flowcellid))
    if df_info.shape[0] == 0 or m3_info.shape[0] == 0 : init("No matching data found.")

    if len(inclusion) > 0:
        print ("inclusion sample:" + "\n".join(inclusion))
        df_info = df_info[ df_info['SAMPLE_ID'].isin(inclusion)]
        if df_info.shape[0] == 0 : init("No corresponding sample IDs.")

    if len(exclusion) > 0:
        print ("exclusion sample:" + ",".join(exclusion))
        df_info = df_info[ ~df_info['SAMPLE_ID'].isin(exclusion)]
        if df_info.shape[0] == 0 : init("No corresponding sample IDs.")

    df_info['PRJ_TYPE'] = df_info['PRJ_TYPE'].str.replace('EWES',"eWES")
    if project_type == "both" :
        df_info = df_info[ df_info['PRJ_TYPE'].isin(['eWES','WTS']) ]
    else :
        df_info = df_info[ df_info['PRJ_TYPE'] == project_type ]
    if df_info.shape[0] == 0 : init("Test type error: no sample ID corresponds.")

    df_info = df_info[~df_info['SAMPLE_ID'].str.contains('_PCE_|_NCE_|_PCT_|_NCT_', regex=True, na=False)]
    if df_info.shape[0] == 0 : init("No clinical specimens match the criteria.")

    df_info = df_info[ df_info['SAMPLE_ID'].isin(m3_info['SAMPLE_ID']) ]
    if df_info.shape[0] == 0 : init("No M3 specimens match the criteria.")

    uniq_info = fcDir_table(df_info, directory)
    if uniq_info.shape[0] == 0: init()

    df_info = pd.merge(df_info, uniq_info, on=['sub_name','PRJ_TYPE'])
    os.makedirs(outdir, exist_ok=True)

    for pj_type in df_info['PRJ_TYPE'].unique():

        df_prj = df_info[df_info['PRJ_TYPE']==pj_type].reset_index(drop=True)
        anal_dir = os.path.join(directory,pj_type,df_prj['seqDir'][0])

        if pj_type == 'eWES':
            out_file_1 = os.path.join(outdir, df_prj['seqDir'][0] + '.SNV_INDEL.tsv')
            out_file_2 = os.path.join(outdir, df_prj['seqDir'][0] + '.CSV.tsv')
            out_file_3 = os.path.join(outdir, df_prj['seqDir'][0] + '.TMB_MSI.tsv')
            remove_files([out_file_1,out_file_2])

            df_SNV = pd.DataFrame(columns=["Sample", "Diagnosis"] + use_columns_ewes_1)
            df_CNV = pd.DataFrame(columns=["Sample", "Diagnosis"] + use_columns_ewes_2)
            df_TMB_MSI = pd.DataFrame(columns=["Sample", "Diagnosis", "TMB score", "TMB status", "MSI score", "MSI status"])

            for i, item in df_prj.iterrows() :

                # SNV & InDel
                file_snv = os.path.join(anal_dir, item['SAMPLE_ID'], 'Summary', item['SAMPLE_ID']+'.summarized.snv.target.tsv')
                if not os.path.isfile(file_snv):
                    print('file not exists: ' + file_snv)
                else:
                    data = pd.read_csv(file_snv, sep="\t")
                    data = data[use_columns_ewes_1].drop_duplicates()
                    data['Clinvar_CLNSIG'] = data['Clinvar_CLNSIG'].astype(str).replace("nan", None)
                    data['ONCOKB_ONCOGENICITY'] = data['ONCOKB_ONCOGENICITY'].astype(str).replace("nan", None)
                    data['REPORT'] = np.where(data['Clinvar_CLNSIG'].str.contains('Pathogenic|Likely_pathogenic', case=True, na=False),'PASS',
                                np.where(data['ONCOKB_ONCOGENICITY'].str.contains('oncogenic|resistance', case=False, na=False),'PASS',''))
                    data = pd.merge(data, var_tbl, how='left', on='SYMBOL')
                    if data.shape[0] == 0 :
                        data.loc[0] = ['-'] * len(use_columns_ewes_1)

                    data.insert(0,'Diagnosis',item['DIAGNOSIS_NAME'])
                    data.insert(0,'Sample',item['SAMPLE_ID'])

                    if df_SNV.empty:
                        df_SNV = data.copy()
                    else:
                        df_SNV = pd.concat([df_SNV, data], axis=0, ignore_index=True)

                # CNV
                file_cnv = os.path.join(anal_dir, item['SAMPLE_ID'], 'CNV', item['SAMPLE_ID']+'.cnv.marked.tsv')
                if not os.path.isfile(file_cnv):
                    print('file not exists: ' + file_cnv)
                else :
                    data = pd.read_csv(file_cnv, sep="\t", low_memory=False)
                    data = data[ data['FILTER_GENES']=='PASS' ]
                    data = data[ data['FILTER_CN']=='PASS' ]
                    data = data[use_columns_ewes_2].drop_duplicates()
                    if data.shape[0] == 0 :
                        data.loc[0] = ['-'] * len(use_columns_ewes_2)

                    data.insert(0, 'Diagnosis', item['DIAGNOSIS_NAME'])
                    data.insert(0, 'Sample', item['SAMPLE_ID'])

                    if df_CNV.empty:
                        df_CNV = data.copy()
                    else:
                        df_CNV = pd.concat([df_CNV, data], axis=0, ignore_index=True)

                # TMB & MSI
                TMB_status = pic_value(os.path.join(anal_dir, item['SAMPLE_ID'], 'Summary', item['SAMPLE_ID']+'.summarized.tmb.exome.tsv'), 'TMB_STATUS')
                MSI_status = pic_value(os.path.join(anal_dir, item['SAMPLE_ID'], 'Summary', item['SAMPLE_ID']+'.summarized.msi.exome.tsv'), 'MSI_STATUS')
                MSI_status = 'MSI-H not detected' if MSI_status == 'MSS' else MSI_status
                TMB_score = pic_value(os.path.join(anal_dir, item['SAMPLE_ID'], 'Summary', item['SAMPLE_ID']+'.summarized.tmb.exome.tsv'), 'TMB')
                MSI_score = pic_value(os.path.join(anal_dir, item['SAMPLE_ID'], 'Summary', item['SAMPLE_ID']+'.summarized.msi.exome.tsv'), 'MSI')

                df_TMB_MSI.loc[len(df_TMB_MSI)] = [item['SAMPLE_ID'],item['DIAGNOSIS_NAME'],TMB_score,TMB_status,MSI_score,MSI_status ]

            df_SNV.to_csv(out_file_1, sep="\t", index=False)
            df_CNV.to_csv(out_file_2, sep="\t", index=False)
            df_TMB_MSI.to_csv(out_file_3, sep="\t", index=False)

        elif pj_type == 'WTS':
            out_file_1 = os.path.join(outdir, df_prj['seqDir'][0] + '.Fusion.tsv')
            out_file_2 = os.path.join(outdir, df_prj['seqDir'][0] + '.Skipped.tsv')
            remove_files([out_file_1,out_file_2])

            df_fusion = pd.DataFrame(columns=["Sample", "Diagnosis"] + use_columns_wts_1)
            df_skipped = pd.DataFrame(columns=["Sample", "Diagnosis"] + use_columns_wts_2)

            for i, item in df_prj.iterrows() :

                # Fusion
                file_mk = os.path.join(anal_dir, item['SAMPLE_ID'], 'Fusion', item['SAMPLE_ID']+'.fusion.marked.tsv')
                if not os.path.isfile(file_mk):
#                    print('file not exists: ' + file_mk)
                    data_f = pd.DataFrame(columns=use_columns_wts_1)
                else:
                    data_f = pd.read_csv(file_mk, sep="\t", low_memory=False)
                    data_f = data_f[['FILTER_ONCOKB','gene1','gene2','chr1','breakpoint_1','chr2','breakpoint_2','max_split_cnt','max_span_cnt','sample_type','disease','tools','inferred_fusion_type','samples','cancer_db_hits','fusion_IDs']]
                    data_f = data_f.drop_duplicates()
                    data_f = expand_breakpoints(data_f)
                    data_f = data_f.rename(columns={'FILTER_ONCOKB':'OncoKB'})
                    data_f.insert(1,'cancer-related',"")
                    data_f.insert(0,'Out-of-Frame',"PASS")
                    data_f['cancer_db_hits'] = data_f['cancer_db_hits'].astype(str)

                file_cis = os.path.join(anal_dir, item['SAMPLE_ID'], 'Fusion', 'Metafusion', 'final.n2.cluster.CANCER_FUSIONS.cis-sage.filtered')
                if not os.path.isfile(file_cis):
                    print('file not exists: ' + file_cis)
                    data = pd.DataFrame(columns=use_columns_wts_1)
                else :
                    data = pd.read_csv(file_cis, sep="\t", low_memory=False)
                    data = data.rename(columns={'#gene1':'gene1'})
                    data['cancer_db_hits'] = data['cancer_db_hits'].astype(str)
                    data = expand_breakpoints(data)

                if data.shape[0] > 0:
                    data = pd.merge(data, data_f, on=data.columns.tolist(), how='outer')[use_columns_wts_1]
                    data['Out-of-Frame'] = data['Out-of-Frame'].fillna('FAIL')
                    data = data.sort_values(['gene1','gene2','breakpoint_1','breakpoint_2'])
                elif data_f.shape[0] > 0 :
                    data = data_f.copy()
                else :
#                    data.loc[0] = ['-'] * len(use_columns_wts_1) 
                    data['samples'] = item['SAMPLE_ID']

                data.insert(1, 'Diagnosis', item['DIAGNOSIS_NAME'])

                if df_fusion.empty:
                    df_fusion = data.copy()
                else :
                    df_fusion = pd.concat([df_fusion, data], axis=0, ignore_index=True)

                # Exon Skipped
                file_es = os.path.join(anal_dir, item['SAMPLE_ID'], 'Alternative_splicing', 'ESDetector', item['SAMPLE_ID']+'.exon_skipped.tsv')
                if not os.path.isfile(file_es):
                    print('file not exists: ' + file_es)
                else:
                    data = pd.read_csv(file_es, sep="\t")
                    if data.shape[0] > 0 :

                        data = data.rename(columns={'#gene1':'gene1'})
                        data.insert(0, 'Diagnosis', item['DIAGNOSIS_NAME'])
                        data.insert(0, 'Sample', item['SAMPLE_ID'])

                        if df_skipped.empty:
                            df_skipped = data.copy()
                        else:
                            df_skipped = pd.concat([df_skipped, data], axis=0, ignore_index=True)

            df_skipped['ratio'] = pd.to_numeric(df_skipped['ratio'], errors="coerce")
            df_skipped['tpm_variant'] = pd.to_numeric(df_skipped['tpm_variant'], errors="coerce")
            df_skipped['Report'] = df_skipped.apply(lambda row: "PASS" if (row['ratio'] >= 1) and (row['tpm_variant'] >= 0.01) else "FAIL", axis=1)
            df_skipped['Review'] = ''

            df_fusion.to_csv(out_file_1, sep="\t", index=False)
            df_skipped.to_csv(out_file_2, sep="\t", index=False)

        else :
            print('PRJ_TYPE error:' + pj_type )
            continue


