import pandas as pd
from .func import *


def load_value(file):
    try :
        data = pd.read_csv(file, header=0)
        return data[['Purity','Ploidy']].values.tolist()[0]
    except Exception as e:
        return [ 'NA', 'NA' ]

def run_cnv(args):

    flowcellid = args.flowcellid
    directory = args.directory
    outdir = args.outdir
    inclusion = [x.strip() for x in args.inclusion.split(',') if not x.strip() == '']
    exclusion = [x.strip() for x in args.exclusion.split(',') if not x.strip() == '']

    inclusion = rmdup_list(inclusion)
    exclusion = rmdup_list(exclusion)

    if len(inclusion) > 0 and len(exclusion) > 0:
        print('ERROR: Inclusion and exclusion cannot be specified simultaneously.')
        parser.print_help()
        sys.exit(1)

    df_info = getinfo(SelectData(flowcellid))
    if df_info.shape[0] == 0 : init()

    if len(inclusion) > 0:
        print ("inclusion sample:" + "\n".join(inclusion))
        df_info = df_info[ df_info['SAMPLE_ID'].isin(inclusion)]
        if df_info.shape[0] == 0 : init("No corresponding sample IDs.")

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
    anal_dir = os.path.join(directory,"eWES",df_info['seqDir'][0])

    out_df = pd.DataFrame(columns=['sample_id','bin_size','purity_400','ploidy_400','purity_800','ploidy_800','purity_1600','ploidy_1600'])
    for i, item in df_info.iterrows() :
        file = os.path.join(anal_dir, item['SAMPLE_ID'], 'CNV', 'PureCN', 'bin_size.txt')
        try :
            f = open(file, 'r')
            fac = [item['SAMPLE_ID'], f.read()]
            f.close()
        except Exception as e:
            print('CNV process has not been completed:' + item['SAMPLE_ID'])
            fac = [item['SAMPLE_ID'], '-']

        for bin in ['400','800','1600']:
            file = os.path.join(anal_dir, item['SAMPLE_ID'], 'CNV', 'PureCN', bin, item['SAMPLE_ID'] + '.tumour.exome.purecn.csv')
            fac.extend( load_value(file) )

        out_df.loc[out_df.shape[0]+1] = fac

    os.makedirs(outdir, exist_ok=True)
    out_df.to_csv(os.path.join(outdir, df_info['seqDir'][0] + '.tsv'), header=True, index=False, sep="\t")


