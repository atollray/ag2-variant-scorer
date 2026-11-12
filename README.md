# ag2-variant-scorer

AlphaGenome 2の予測APIで公開変異を評価し、同じ変異のスコアを細胞条件ごとに記録するためのローカル検証プロジェクト。

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

## 再利用したコード

`reference.py` のFASTA読み取りとSNVのREF照合は、`atollray/alphagenome-variant-scorer` のコミット `9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc` にある `Reference` と `verified_normalize` の該当部分を基にしている。同コミットのApache-2.0ライセンスを `LICENSE` に収録した。[FASTA読み取りのコピー元](https://github.com/atollray/alphagenome-variant-scorer/blob/9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc/src/alphagenome_variant_scorer/reference.py)・[REF照合のコピー元](https://github.com/atollray/alphagenome-variant-scorer/blob/9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc/src/alphagenome_variant_scorer/variants.py)
