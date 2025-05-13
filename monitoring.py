import argparse
from modules import *

VERSION="v1.0.0"

class CustomHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        version_info = f"version: {VERSION}\n"
        self._add_item(lambda: version_info, [])
        super().add_usage(usage, actions, groups, prefix)

def main():
    parser = argparse.ArgumentParser(
        description="Tools for monitoring analysis data.",
        formatter_class=CustomHelpFormatter
    )
    parser.add_argument('--version','-v', action='version', version=f'%(prog)s {VERSION}')
    subparsers = parser.add_subparsers(dest="command", required=True)

    # monitoring QC
    parser_qc = subparsers.add_parser("QC", help="QC monitoring")
    parser_qc.add_argument("--output","-o", required=False, help="output file path", default="")
    parser_qc.set_defaults(func=run_qc)

    # monitoring CNV (pureCN)
    parser_cnv = subparsers.add_parser("CNV", help="CNV(PureCN) monitoring", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_cnv.add_argument("--flowcellid","-fc", required=True, help="flowcell id")
    parser_cnv.add_argument("--inclusion","-i", required=False, help="sample IDs to include (comma separated)", default="")
    parser_cnv.add_argument("--exclusion","-e", required=False, help="sample IDs to exclude (comma separated)", default="")
    parser_cnv.add_argument("--directory","-d", required=False, help="parent analytical directory", default="/data1/data/result")
    parser_cnv.add_argument("--outdir","-o", required=False, help="output directory path", default="/data1/work/monitoring/PureCN")
    parser_cnv.set_defaults(func=run_cnv)

    # monitoring Fusion (STAR-SEQR)
    parser_seqr = subparsers.add_parser("fusion", aliases=['FS'], help="Fusion(STAR-SEQR) monitoring", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_seqr.add_argument("--sample","-s", required=True, help="sample id")
    parser_seqr.add_argument("--verbose","-v", required=False, help="Show details", action='store_true')
    parser_seqr.add_argument("--analysis_dir","-d", required=False, help="parent analytical directory", default="/data1/data/result")
    parser_seqr.set_defaults(func=run_seqr)

    # create pre-filtered data
    parser_pre = subparsers.add_parser("preFilter", aliases=['PRE'], help="create pre-filtered data", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_pre.add_argument("--flowcellid","-fc", required=True, help="flowcell id")
    parser_pre.add_argument("--directory","-d", required=False, help="parent analytical directory", default="/data1/data/result")
    parser_pre.add_argument("--project_type","-t", required=False, help="project type", default="both", choices=["both","WTS","eWES"])
    parser_pre.add_argument("--outdir","-o", required=False, help="output directory path", default="/data1/work/monitoring/preFilter")
    parser_pre.add_argument("--inclusion","-i", required=False, help="sample IDs to include (comma separated)", default="")
    parser_pre.add_argument("--exclusion","-e", required=False, help="sample IDs to exclude (comma separated)", default="")
    parser_pre.set_defaults(func=run_preFilter)

    # Creating Benchmark Files
    parser_bcm = subparsers.add_parser("benchmark", aliases=['BM'], help="List benchmark data.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_bcm.add_argument("--flowcellid","-fc", required=True, help="flowcell id")
    parser_bcm.add_argument("--project_type","-t", required=False, help="project type", default="both", choices=["both","WTS","eWES"])
    parser_bcm.add_argument("--directory","-d", required=False, help="parent analytical directory", default="/data1/data/result")
    parser_bcm.add_argument("--outdir","-o", required=False, help="output directory path", default="/data1/work/monitoring/benchmark")
    parser_bcm.add_argument("--inclusion","-i", required=False, help="sample IDs to include (comma separated)", default="")
    parser_bcm.add_argument("--exclusion","-e", required=False, help="sample IDs to exclude (comma separated)", default="")
    parser_bcm.set_defaults(func=run_benchmark)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":

    main()

