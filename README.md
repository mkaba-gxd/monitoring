# モニタリング
解析済（または解析中）データのモニタリングを行う。

| command        | 概要                             |
|:---------------|:---------------------------------|
| QC             | WET,DRY工程のQC値一覧作成         |
| CNV            | (PureCN) purity, ploidyの一覧作成 |
| fusion, FS     | (STAR-SEQR) 所要時間の推定        |
| splice, AS     | EGFR, MET, AR 領域のdepthを描画   |
| preFilter, PRE | フィルター前データ作成             |
| benchmark, BM  | 工程所要時間の一覧作成             |

## 変数の定義(共通)
```bash
img=/data1/labTools/labTools.sif
SCRIPT=/data1/labTools/monitoring/latest/monitoring.py
```
## マニュアルの表示
全体の概要表示
```bash
$ singularity exec --bind /data1 $img python $SCRIPT --help
version: v1.0.0
usage: monitoring.py [-h] [--version] {QC,CNV,fusion,FS,splice,AS,preFilter,PRE,benchmark,BM} ...

Tools for monitoring analysis data.

positional arguments:
  {QC,CNV,fusion,FS,splice,AS,preFilter,PRE,benchmark,BM}
    QC                  QC monitoring
    CNV                 CNV(PureCN) monitoring
    fusion (FS)         Fusion(STAR-SEQR) monitoring
    splice (AS)         Alternative Splicing monitoring
    preFilter (PRE)     create pre-filtered data
    benchmark (BM)      List benchmark data.

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit
```
コマンド別の詳細表示
```
singularity exec --bind /data1 $img python $SCRIPT <command> --help
```

## 1\. QC値
スクリプトを実行した時点でデータベースに登録されている全検体のQC情報の一覧を作成する。\
QC項目はOncoStationで確認できる項目と同じ。
```
singularity exec --bind /data1 $img python $SCRIPT QC
```
⇒ /data1/work/monitoring/QC/[timestamp].xlsx が作成される

## 2\. CNV（Copy Number Variants/PureCN）
PC,NCを除いた各サンプルについて、解析で採用されたbin sizeと bin size 400/800/1600 のPureCNで算出されたpurityとploidyの一覧を出力する。
### オプションの詳細
```
$ singularity exec --bind /data1 $img python $SCRIPT CNV --help
version: v1.0.0
usage: monitoring.py CNV [-h] --flowcellid FLOWCELLID [--inclusion INCLUSION] [--exclusion EXCLUSION]
                         [--directory DIRECTORY] [--outdir OUTDIR]
optional arguments:
  -h, --help            show this help message and exit
  --flowcellid FLOWCELLID, -fc FLOWCELLID
                        flowcell id (default: None)
  --inclusion INCLUSION, -i INCLUSION
                        sample IDs to include (comma separated) (default: )
  --exclusion EXCLUSION, -e EXCLUSION
                        sample IDs to exclude (comma separated) (default: )
  --directory DIRECTORY, -d DIRECTORY
                        parent analytical directory (default: /data1/data/result)
  --outdir OUTDIR, -o OUTDIR
                        output directory path (default: /data1/work/monitoring/PureCN)
```
| option          |required | 概要                     |default                       |
|:----------------|:-------:|:-------------------------|:-----------------------------|
|--flowcellid/-fc |True     |バッチ固有のID。OncoStationに掲載されている9桁の半角英数字   |None |
|--inclusion/-i   |False    |出力するSample IDを限定。カンマ区切りで複数指定可能 |None    |
|--exclusion/-e   |False    |除外するSample IDを指定。カンマ区切りで複数指定可能 |None    |
|--directory/-d   |False    |解析フォルダの親ディレクトリ |/data1/data/result            |
|--outdir/-o      |False    |結果の出力先ディレクトリ     |/data1/work/monitoring/PureCN |

⇒ /data1/work/monitoring/PureCN/[batchfolder].tsv が作成される。\
*すでに出力ファイルが存在する場合は上書きする。

## 3\. fusion（STAR-SEQR）
STAR-RSEQの実行時間の目安となる sequenceの組合せ総数を算出する。\
値が 10^6 未満なら数時間で終了する可能性が高い。
### オプションの詳細
```
$ singularity exec --bind /data1 $img python $SCRIPT fusion --help
version: v1.0.0
usage: monitoring.py fusion [-h] --sample SAMPLE [--verbose] [--analysis_dir ANALYSIS_DIR]

optional arguments:
  -h, --help            show this help message and exit
  --sample SAMPLE, -s SAMPLE
                        sample id (default: None)
  --verbose, -v         Show details (default: False)
  --analysis_dir ANALYSIS_DIR, -d ANALYSIS_DIR
                        parent analytical directory (default: /data1/data/result)
```
| option           |required | 概要                     |default             |
|:-----------------|:-------:|:-------------------------|:-------------------|
|--sample/-s       |True     |Sample ID。複数指定不可    |None                |
|--verbose/-v      |False    |詳細を表示するかどうか      |False               |
|--analysis_dir/-d |False    |解析フォルダの親ディレクトリ |/data1/data/result |

⇒ sequenceの組合せ総数がディスプレイに表示される\
2025/4/24 時点：組合せ総数の最大11,049,185に対し、STAR-SEQRの所要時間は 44:29:04\
2025/6/2 時点：組合せ総数の最大30,473,853に対し、STAR-SEQRの所要時間は 166:37:57

## 4\. splice（Alternative Splicing）
BAMファイルからEGFR, MET,AR領域のdepthを計測し、exon領域とともに描画する。
### オプションの詳細
```
$ singularity exec --bind /data1 $img python $SCRIPT splice --help
version: v1.0.0
usage: monitoring.py splice [-h] --sample SAMPLE [--category CATEGORY]
                            [--analysis_dir ANALYSIS_DIR] [--outdir OUTDIR]
optional arguments:
  -h, --help            show this help message and exit
  --sample SAMPLE, -s SAMPLE
                        sample id (default: None)
  --category CATEGORY, -c CATEGORY
                        Comma-separated list of genes to include. Valid options: EGFR, MET, AR (default: ['AR'])
  --analysis_dir ANALYSIS_DIR, -d ANALYSIS_DIR
                        parent analytical directory (default: /data1/data/result)
  --outdir OUTDIR, -o OUTDIR
                        output directory path (default: /data1/work/monitoring/splice)
```
| option           |required | 概要                     |default                   |
|:-----------------|:-------:|:-------------------------|:-------------------------|
|--sample/-s       |True     |Sample ID, 複数指定不可    |None                      |
|--category/c      |False    |描画する遺伝子領域。EGFR, MET, ARから選択 | AR          |
|--analysis_dir/-d |False    |解析フォルダの親ディレクトリ |/data1/data/result        |
|--outdir/-o       |False    |結果の出力先ディレクトリ |/data1/work/monitoring/splice |

⇒ /data1/work/monitoring/splice/[sample]\_dnacopy\_[EGFR/MET/AR].pdf が作成される。

## 5\. preFilter
Filter前の解析データをExcel出力する。
### オプションの詳細
```
$ singularity exec --bind /data1 $img python $SCRIPT preFilter --help
version: v1.0.0
usage: monitoring.py preFilter [-h] --flowcellid FLOWCELLID [--directory DIRECTORY] [--project_type {both,WTS,eWES}]
                               [--outdir OUTDIR] [--inclusion INCLUSION] [--exclusion EXCLUSION]
optional arguments:
  -h, --help            show this help message and exit
  --flowcellid FLOWCELLID, -fc FLOWCELLID
                        flowcell id (default: None)
  --directory DIRECTORY, -d DIRECTORY
                        parent analytical directory (default: /data1/data/result)
  --project_type {both,WTS,eWES}, -t {both,WTS,eWES}
                        project type (default: both)
  --outdir OUTDIR, -o OUTDIR
                        output directory path (default: /data1/work/monitoring/preFilter)
  --inclusion INCLUSION, -i INCLUSION
                        sample IDs to include (comma separated) (default: )
  --exclusion EXCLUSION, -e EXCLUSION
                        sample IDs to exclude (comma separated) (default: )
```
| option           |required | 概要                     |default                           |
|:-----------------|:-------:|:-------------------------|:---------------------------------|
|--flowcellid/-fc  |True     |バッチ固有のID。OncoStationに掲載されている9桁の半角英数字 |None |
|--directory/-d    |False    |解析フォルダの親ディレクトリ |/data1/data/result                |
|--project_type/-t |False    |解析種別。both,eWES,WTSから選択                   |both        |
|--outdir/-o       |False    |結果の出力先ディレクトリ     |/data1/work/monitoring/preFilter |
|--inclusion/-i    |False    |出力するSample IDを限定。カンマ区切りで複数指定可能 |None        |
|--exclusion/-e    |False    |除外するSample IDを指定。カンマ区切りで複数指定可能 |None        |

⇒ /data1/work/monitoring/preFilter/[batchfolder] の下に複数の.xlsxファイルが作成される

## 6\. benchmark
解析工程でBenchmarkフォルダに出力される各工程の所要時間(h:m:sの値)のテーブルをファイル出力する。
### オプションの詳細
```
$ singularity exec --bind /data1 $img python $SCRIPT benchmark --help
version: v1.0.0
usage: monitoring.py benchmark [-h] --flowcellid FLOWCELLID [--project_type {both,WTS,eWES}] [--directory DIRECTORY]
                               [--outdir OUTDIR] [--inclusion INCLUSION] [--exclusion EXCLUSION]
optional arguments:
  -h, --help            show this help message and exit
  --flowcellid FLOWCELLID, -fc FLOWCELLID
                        flowcell id (default: None)
  --project_type {both,WTS,eWES}, -t {both,WTS,eWES}
                        project type (default: both)
  --directory DIRECTORY, -d DIRECTORY
                        parent analytical directory (default: /data1/data/result)
  --outdir OUTDIR, -o OUTDIR
                        output directory path (default: /data1/work/monitoring/benchmark)
  --inclusion INCLUSION, -i INCLUSION
                        sample IDs to include (comma separated) (default: )
  --exclusion EXCLUSION, -e EXCLUSION
                        sample IDs to exclude (comma separated) (default: )
```
| option           |required | 概要                     |default                           |
|:-----------------|:-------:|:-------------------------|:---------------------------------|
|--flowcellid/-fc  |True     |バッチ固有のID。OncoStationに掲載されている9桁の半角英数字 |None |
|--project_type/-t |False    |解析種別。both,eWES,WTSから選択                   |both        |
|--directory/-d    |False    |解析フォルダの親ディレクトリ |/data1/data/result                |
|--outdir/-o       |False    |結果の出力先ディレクトリ     |/data1/work/monitoring/benchmark |
|--inclusion/-i    |False    |出力するSample IDを限定。カンマ区切りで複数指定可能 |None        |
|--exclusion/-e    |False    |除外するSample IDを指定。カンマ区切りで複数指定可能 |None        |

⇒ /data1/work/monitoring/benchmark/[batchfolder].xlsx が作成される\
*すでに出力ファイルが存在する場合は上書きする。
