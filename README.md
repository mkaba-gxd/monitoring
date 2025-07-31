# monitoring ツール
CAP検査（eWES/WTS）で実施された解析について、モニタリングを行う。\
指定されたsample IDやflowcell IDから検体情報をデータベースに問合せ、解析ディレクトリ内のデータをロードするため、データベースに登録がない検体や、規程の場所に解析データやログファイルがない検体に対しては実行できません。\
**また、データベースの設計内容が不明なため、データベース検索時に想定外の動作を行う可能性があります。**
| command           | 概要                              |
|:------------------|:----------------------------------|
| [QC](#QC)         | WET,DRY工程のQC値一覧作成          |
| [pureCN](#pureCN) | PureCNで算出されたpurity, ploidyの一覧作成 |
| [CNV](#CNV)       | 指定した遺伝子セットのコピー数一覧を作成 |
| [SNV](#SNV)       | 指定した場所で観測された変異を抽出 |
| [fusion, FS](#fusion)        | (STAR-SEQR) 所要時間の推定         |
| [splice, AS](#splice)        | EGFR, MET, AR 領域のdepthを描画    |
| [preFilter, PRE](#preFilter) | フィルター前データ作成              |
| [benchmark, BM](#benchmark)  | 工程所要時間の一覧作成              |

## エイリアスの作成 ※初回のみ
~/bin フォルダ直下に以下のコマンドを記載したテキストファイル monitoring を作成し、実行権限を付与する。\
(gxd_pipeline, guest_user ユーザーには実装済み) \
エイリアスを作成しない場合は、singularity でコンテナとスクリプトファイルを指定して実行する。
```
singularity exec --disable-cache --bind /data1 /data1/labTools/labTools.sif python /data1/labTools/monitoring/latest/monitoring.py $@
```
helpページを表示してエイリアスの設定を確認する。以下が表示されればOK。
```
$ monitoring --help
version: v1.1.0
usage: monitoring.py [-h] [--version] {QC,pureCN,CNV,SNV,fusion,FS,splice,AS,preFilter,PRE,benchmark,BM} ...

Tools for monitoring analysis data.

positional arguments:
  {QC,pureCN,CNV,SNV,fusion,FS,splice,AS,preFilter,PRE,benchmark,BM}
    QC                  QC monitoring
    pureCN              PureCN window size monitoring
    CNV                 Copy Numver monitoring
    SNV                 Extract SNV intermediate data.
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
monitoring <command> --help
```

<a id="QC"></a>
## 1\. QC値
スクリプトを実行した時点でデータベースに登録されている全検体のQC情報の一覧を作成する。\
QC項目はOncoStationで確認できる項目と同じ。
```
monitoring QC
```
⇒ /data1/work/monitoring/QC/[timestamp].xlsx が作成される

<a id="pureCN"></a>
## 2\. pureCN
PC,NCを除いた各サンプルについて、解析で採用されたbin sizeと、PureCNで算出された bin size 400/800/1600 のpurityとploidyの一覧を出力する。
```
monitoring pureCN --flowcellid <flowcellid>
```
⇒ /data1/work/monitoring/PureCN/[batchfolder].tsv が作成される。\
※ すでに出力ファイルが存在する場合は上書きする。
<details>
  <summary> 
    More Details
  </summary>

### オプションの詳細
```
$ monitoring pureCN --help
version: v1.1.0
usage: monitoring.py pureCN [-h] --flowcellid FLOWCELLID [--inclusion INCLUSION]
                [--exclusion EXCLUSION] [--directory DIRECTORY] [--outdir OUTDIR]

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
|--flowcellid/-fc |True     |バッチ固有のID。OncoStationに掲載されている9桁の半角英数字 |None |
|--inclusion/-i   |False    |出力するSample IDを限定。カンマ区切りで複数指定可能 |None  |
|--exclusion/-e   |False    |除外するSample IDを指定。カンマ区切りで複数指定可能 |None  |
|--directory/-d   |False    |解析フォルダの親ディレクトリ |/data1/data/result           |
|--outdir/-o      |False    |結果の出力先ディレクトリ     |/data1/work/monitoring/PureCN |

</details>

<a id="CNV"></a>
## 3\. CNV
指定した遺伝子セットのコピー数一覧を出力する。\
※ スクリプトを実行した時点で解析が終了しているすべてのM3検体(Patient IDがM3から始まるもの)のコピー数を検索する。\
※ 解析フォルダに保存されている中間ファイルが削除されている場合は除外される。
```
monitoring CNV --genes [GENES/gene list filepath]
```
⇒ /data1/work/monitoring/CNV/[timestamp].xlsx が作成される。
<details>
  <summary> 
    More Details
  </summary>

### オプションの詳細
```
$ monitoring CNV --help
version: v1.1.0
usage: monitoring.py CNV [-h] --genes GENES [--exclusion EXCLUSION] [--directory DIRECTORY] [--outdir OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  --genes GENES, -g GENES
                        gene names (comma separated) or file paths for a list of gene sets (default: )
  --exclusion EXCLUSION, -e EXCLUSION
                        sample IDs to exclude (comma separated) (default: )
  --directory DIRECTORY, -d DIRECTORY
                        parent analytical directory (default: /data1/data/result)
  --outdir OUTDIR, -o OUTDIR
                        output directory path (default: /data1/work/monitoring/CNV)
```
| option        |required | 概要                     |default                     |
|:--------------|:-------:|:-------------------------|:--------------------------|
|--genes/-g     |True     |遺伝子名(カンマ区切りで複数指定可能)または遺伝子名を記載したファイルパス |None |
|--exclusion/-e |False    |除外するSample IDを指定。カンマ区切りで複数指定可能 |None |
|--directory/-d |False    |解析フォルダの親ディレクトリ |/data1/data/result         |
|--outdir/-o    |False    |結果の出力先ディレクトリ     |/data1/work/monitoring/CNV |

</details>

<a id="SNV"></a>
## 4. SNV
```
monitoring SNV --sample <sample> --position <chr:pos> (--window <int>)
```
⇒ 指定した場所について、mutect2,lofreq,freebayesで検出された変異を表示する。

<details>
  <summary>
    More Details
  </summary>

### オプションの詳細
```
$ monitoring SNV --help
version: v1.1.0
usage: monitoring.py SNV [-h] --sample SAMPLE --position POSITION [--window WINDOW] [--directory DIRECTORY]

optional arguments:
  -h, --help            show this help message and exit
  --sample SAMPLE, -s SAMPLE
                        sample id (default: None)
  --position POSITION, -p POSITION
                        site of mutation. (default: None)
  --window WINDOW, -w WINDOW
                        locus width. (default: 0)
  --directory DIRECTORY, -d DIRECTORY
                        parent analytical directory (default: /data1/data/result)

```
| option        |required | 概要                     |default               |
|:--------------|:-------:|:-------------------------|:---------------------|
|--sample/-s    |True     |Sample ID。複数指定不可   |None                  |
|--position/-v  |True     |検索するゲノム上の場所(chr:position) |None       |
|--window/-w    |False    |検索する範囲(positionの前後<window>bpを検索) |0  |
|--directory/-d |False    |解析フォルダの親ディレクトリ |/data1/data/result |

</details>

<a id="fusion"></a>
## 5\. fusion（STAR-SEQR）
STAR-RSEQの実行時間の目安となる sequenceの組合せ総数を算出する。\
値が 10^6 未満なら数時間で終了する可能性が高い。
```
monitoring fusion --sample <sample>
monitoring FS -s <sample>
```
⇒ sequenceの組合せ総数がディスプレイに表示される\
2025/4/24 時点：組合せ総数の最大11,049,185に対し、STAR-SEQRの所要時間は 44:29:04\
2025/6/2 時点：組合せ総数の最大30,473,853に対し、STAR-SEQRの所要時間は 166:37:57
<details>
  <summary> 
    More Details
  </summary>
  
### オプションの詳細
```
$ monitoring fusion --help
version: v1.1.0
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

</details>

<a id="splice"></a>
## 6\. splice（Alternative Splicing）
BAMファイルからEGFR, MET,AR領域のdepthを計測し、exon領域とともに描画する。
```
monitoring splice --sample <sample>
monitoring AS -s <sample>
```
⇒ /data1/work/monitoring/splice/[sample]\_dnacopy\_[EGFR/MET/AR].pdf が作成される。
<details>
  <summary> 
    More Details
  </summary>

### オプションの詳細
```
$ monitoring splice --help
version: v1.1.0
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

</details>

<a id="preFilter"></a>
## 7\. preFilter
Filter前の解析データをExcel出力する。
```
monitoring preFilter --flowcellid <flowcellid>
monitoring PRE -fc <flowcellid>
```
⇒ /data1/work/monitoring/preFilter/[batchfolder] の下に複数の.xlsxファイルが作成される
<details>
  <summary> 
    More Details
  </summary>

### オプションの詳細
```
$ monitoring preFilter --help
version: v1.1.0
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

</details>

<a id="benchmark"></a>
## 8\. benchmark
解析工程でBenchmarkフォルダに出力される各工程の所要時間(h:m:sの値)のテーブルをファイル出力する。
```
monitoring benchmark --flowcellid <flowcellid>
monitoring BM -fc <flowcellid>
```
⇒ /data1/work/monitoring/benchmark/[batchfolder].xlsx が作成される\
※ すでに出力ファイルが存在する場合は上書きする。
<details>
  <summary> 
    More Details
  </summary>
  
### オプションの詳細
```
$ monitoring benchmark --help
version: v1.1.0
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

</details>
