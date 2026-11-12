# ag2-variant-scorer

AlphaGenome 2の予測APIで公開SNVと合成構造変異を評価し、変異スコアを細胞条件ごとに記録するためのローカル検証プロジェクト。

## 初回検証の再実行

Python 3.12で、AlphaGenome APIキーと既存の `.fai` を持つGRCh38 FASTAを用意する。参照配列のビルドは出典で確認する。スクリプトは対象座標のREF 1塩基を照合する。PowerShellでリポジトリのルートから実行する。

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe pilot.py --reference 'C:\path\to\hg38.fa' `
  --api-key-file 'C:\path\to\private\api_key.txt' `
  --output-dir 'C:\path\to\data\ag2-pilot-run'
```

キーを `ALPHAGENOME_API_KEY` 環境変数で設定している場合は `--api-key-file` を省略できる。キーファイルと出力先はリポジトリの外に置き、出力先は実行ごとに新しいディレクトリを指定する。

[公式チュートリアル](https://www.alphagenomedocs.com/colabs/v2/conditioning_on_single_cell_data.html)の公開SNV `chr2:219243701:C>T`（`GLB1L`）を使う。REFをローカルで照合し、血液のCD14陽性単球 `UBERON:0000178_CL:0001054` と内皮細胞 `CL:0000115` の2条件で、1,048,576 bpの窓を使って `RNA_SEQ` スコアを取得する。トラックごとの `raw_score` と `imputation_mode` をそのまま記録する。[公式の条件選択ガイド](https://www.alphagenomedocs.com/colabs/v2/choosing_a_condition.html)・[スコアの定義](https://www.alphagenomedocs.com/variant_scoring.html)

| 出力 | 内容 |
| --- | --- |
| `conditions.tsv` | 条件ID、名前、RNA-seqの出力区分 |
| `scores.tsv` | 2条件の遺伝子・トラック別スコア。`requested_condition` で要求条件を識別 |
| `request.json` | 入力、窓、モデル指定、SDK版、取得時刻。両条件の保存後に作成 |

## 合成SVの検証を再実行する

上と同じ実行コマンドに `--sv` を加え、新しい出力先を指定する。GRCh38 `chr22:36202000–36202999`（1-based、両端を含む）の1,000 bp単純欠失と、同方向に1コピーを追加するタンデム重複を評価する。FASTAから両例のREF・ALT配列を構成し、左正規化後にREFと変異後配列を照合する。両例を同じ1,048,576 bp窓と上記2条件でスコアリングする。

`scores.tsv` の `variant_case` で合成例を区別する。`request.json` に、対象区間、正規化後の位置、REF・ALTの全配列、長さ変化を記録する。確認対象はこの2例の `RNA_SEQ` スコア取得。スコアの計算では欠失・挿入後の信号をREF座標に整列する。[公式のスコア定義](https://www.alphagenomedocs.com/variant_scoring.html)

## Splicingの検証を再実行する

同じ実行コマンドに `--splicing --sv` を加え、新しい出力先を指定する。公式のAPOL4 SNV `chr22:36201698:A>C` と上記の合成DEL・タンデムDUPを、S状結腸 `UBERON:0001159` で評価する。`--splicing` のみではAPOL4 SNVを実行する。[公式のAPOL4例](https://www.alphagenomedocs.com/colabs/v2/splicing_variant.html)

| 出力 | 内容 |
| --- | --- |
| 各変異の `scores.tsv` | `SPLICE_SITES`、`SPLICE_SITE_USAGE`、`SPLICE_JUNCTIONS` の `raw_score`。接合部の座標と鎖を保持 |
| 各変異の `REF/` と `ALT/` | 3出力の `values.npz`、列順に対応するトラック表、接合部配列の行順に対応する `junctions.tsv` |
| `summary.tsv` | 変異ごとに、`observed` トラックの成分別最大絶対値と `sites + usage + junctions / 5`。成分が欠ける場合は統合値を空欄とし、欠けた成分を明記 |
| `request.json` | 入力アレル、予測窓、条件、SDK版、REF／ALTそれぞれの区間・解像度・配列形状 |

REF／ALTの位置別予測は、SDKから返った座標でアレル別に保存する。SVの欠失区間・複製区間を比較するときは、構造に応じた座標対応を確かめる。統合値は予測変化の大きさをまとめた値で、各成分の最大値が異なる遺伝子や部位に由来することがある。[公式のスコア定義](https://www.alphagenomedocs.com/variant_scoring.html)

## 再利用したコード

`reference.py` のFASTA読み取りとSNVのREF照合は、`atollray/alphagenome-variant-scorer` のコミット `9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc` にある `Reference` と `verified_normalize` の該当部分を基にしている。同コミットのApache-2.0ライセンスを `LICENSE` に収録した。[FASTA読み取りのコピー元](https://github.com/atollray/alphagenome-variant-scorer/blob/9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc/src/alphagenome_variant_scorer/reference.py)・[REF照合のコピー元](https://github.com/atollray/alphagenome-variant-scorer/blob/9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc/src/alphagenome_variant_scorer/variants.py)
