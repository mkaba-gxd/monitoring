# モニタリング
解析済（または解析中）データのモニタリングを行う。

| command        | 概要                             |
|:---------------|:---------------------------------|
| QC             | WET,DRY工程のQC値一覧作成         |
| CNV            | (PureCN) purity, ploidyの一覧作成 |
| fusion, FS     | (STAR-SEQR) 所要時間の推定        |
| splice, AS     | MET, AR 領域のdepthを描画         |
| preFilter, PRE | フィルター前データ作成             |
| benchmark, BM  | 工程所要時間の一覧作成             |

## 変数の定義(共通)
```bash
img=/data1/labTools/labTools.sif
SCRIPT=/data1/labTools/supplement/monitoring.py
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

## 2\. PureCN
PC,NCを除いた各サンプルについて、解析で採用されたbin sizeと bin size 400/800/1600 のPureCNで算出されたpurityとploidyの一覧を出力する。
変数の設定
```
flowcellid=""
inc_ID=""
exc_ID=""
```
PC,NCを除く全サンプルの情報を書き出す場合
```
singularity exec --bind /data1 $img python $SCRIPT CNV -fc $flowcellid
```
除外するサンプルIDがある場合（comma区切りで複数指定可）
```
singularity exec --bind /data1 $img python $SCRIPT CNV -fc $flowcellid --exclusion $exc_ID
```
限定したいサンプルIDがある場合（comma区切りで複数指定可）
```
singularity exec --bind /data1 $img python $SCRIPT CNV -fc $flowcellid --inclusion $inc_ID
```
⇒ /data1/work/monitoring/PureCN/[batchfolder].tsv が作成される

## 3\. Fusion
STAR-RSEQの実行時間の目安となる sequenceの組合せ総数を算出する。\
値が 10^6 未満なら数時間で終了する可能性が高い。\
2025/4/24 時点での組合せ総数の最大11,049,185に対し、STAR-SEQRの所要時間は 44:29:04\
2025/6/2 時点での組合せ総数の最大30,473,853に対し、STAR-SEQRの所要時間は 166:37:57\
変数の設定
```
sample=""
```
スクリプト実行
```
singularity exec --bind /data1 $img python $SCRIPT FS -s $sample
```
⇒ sequenceの組合せ総数がディスプレイに表示される

## 4\. Alternative Splicing
BAMファイルからMET,AR領域のdepthを計測し、exon領域とともに描画する。\
変数の設定
```
sample=""
```
スクリプト実行
```
singularity exec --bind /data1 $img python $SCRIPT AS -s $sample
```
⇒ /data1/work/monitoring/splice/[sample]\_dnacopy\_[AR/MET].pdf が作成される。

## 5\. pre-Filter
Filer前の解析結果データを作成する。\
変数の設定
```
flowcellid=""
inc_ID=""
exc_ID=""
```
PC,NCを除く全サンプルの情報を書き出す場合
```
singularity exec --bind /data1 $img python $SCRIPT PRE -fc $flowcellid
```
除外するサンプルIDがある場合（comma区切りで複数指定可）
```
singularity exec --bind /data1 $img python $SCRIPT PRE -fc $flowcellid --exclusion $exc_ID
```
限定したいサンプルIDがある場合（comma区切りで複数指定可）
```
singularity exec --bind /data1 $img python $SCRIPT PRE -fc $flowcellid --inclusion $inc_ID
```
⇒ /data1/work/monitoring/preFilter/[batchfolder] の下に複数のxlsxファイルが作成される

## 6\. benchmark
解析工程でBenchmarkフォルダに出力される各工程の所要時間(h:m:sの値)のテーブルを作成する。\
変数の設定
```
flowcellid=""
inc_ID=""
exc_ID=""
```
PC,NCを除く全サンプルの情報を書き出す場合
```
singularity exec --bind /data1 $img python $SCRIPT BM -fc $flowcellid
```
除外するサンプルIDがある場合（comma区切りで複数指定可）
```
singularity exec --bind /data1 $img python $SCRIPT BM -fc $flowcellid --exclusion $exc_ID
```
限定したいサンプルIDがある場合（comma区切りで複数指定可）
```
singularity exec --bind /data1 $img python $SCRIPT BM -fc $flowcellid --inclusion $inc_ID
```
⇒ /data1/work/monitoring/benchmark/[batchfolder].xlsx が作成される
