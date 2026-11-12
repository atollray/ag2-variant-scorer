# ag2-variant-scorer

AlphaGenome 2の予測APIでGRCh38の変異を組織・細胞条件ごとにスコアリングするCLI。取得結果を保存し、中断した実行の再開、別実行でのキャッシュ再利用、発現スコアからの遺伝子候補の集約ができる。

## インストール

Python 3.12以上を使う。取得したソースのルートで、使用するOSに応じてインストールする。

**Windows (PowerShell)**

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install .
.\.venv\Scripts\ag2-variant-scorer.exe --help
```

**Linux**

ディストリビューションが `venv` と `ensurepip` を別パッケージで提供する場合は、使用するPythonに対応する仮想環境用パッケージを用意する。

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install .
.venv/bin/ag2-variant-scorer --help
```

仮想環境を有効にすると `ag2-variant-scorer` で実行できる。以下の例ではこのコマンド名を使う。仮想環境を有効にせず使う場合は、上のOS別の実行パスに置き換える。各サブコマンドの引数は `score --help`、`resume --help`、`summarize --help` で確認できる。

リポジトリからの直接実行にも対応する。従来の `python score.py`、`python score.py --resume DIR`、`python summarize.py` を使う場合は、先に `python -m pip install -r requirements.txt` を実行する。

## 公開例で1件取得する

GRCh38の非圧縮FASTAと、同じ場所にある `<FASTA名>.fai` を用意する。`samtools faidx /path/to/GRCh38.fa` でインデックスを作成できる。FASTAの染色体名は入力の染色体名と揃える。以下のコマンドでは参照ファイルとキーファイルのパスを手元のものに置き換える。

[examples/apol4.tsv](examples/apol4.tsv) は公式チュートリアルの公開APOL4 SNVを使う。S状結腸の `RNA_SEQ` を指定し、まず事前確認を実行する。正しい参照ではREF照合1件、除外0件、予定1件と表示される。この手順はAPIキーなしで実行できる。[公式のAPOL4例](https://www.alphagenomedocs.com/colabs/v2/splicing_variant.html)

```text
ag2-variant-scorer score --input examples/apol4.tsv --reference "/path/to/GRCh38.fa" --condition UBERON:0001159 --scorer RNA_SEQ --dry-run
```

実際に取得するにはAlphaGenome APIキーを用意する。以下の例では出力先に新しい `runs/apol4` を指定する。認証を `ALPHAGENOME_API_KEY` 環境変数で設定済みの場合は `--api-key-file` を省略できる。[AlphaGenome APIキーの取得](https://alphagenome.google/api)

```text
ag2-variant-scorer score --input examples/apol4.tsv --reference "/path/to/GRCh38.fa" --condition UBERON:0001159 --scorer RNA_SEQ --output-dir runs/apol4 --cache-dir cache/ag2 --api-key-file "/path/to/private/api_key.txt"
ag2-variant-scorer summarize --run-dir runs/apol4 --output-dir runs/apol4-summary
```

## 入力形式と複数条件での取得

入力TSVはUTF-8のタブ区切りで、`input_id`、`chromosome`、`position`、`ref`、`alt` を用意する。`input_id` は入力内で一意にする。GRCh38の1-based座標と、REF・ALTの非空のA/C/G/T配列を指定する。indelやSVは、VCF形式のアンカー塩基を含む実配列で表す。スコア取得前にREFの全塩基をFASTAと照合する。

```tsv
input_id	chromosome	position	ref	alt
apol4	chr22	36201698	A	C
```

次のPowerShell例では、黄斑の錐体と網膜色素上皮（RPE）について、細胞タイプ別の遺伝子発現、バルクRNA-seq、スプライシング3種を取得する。`--condition` にはAG2の参照条件CURIE、`--scorer` には取得するスコア名を1項目ずつ指定する。同じスコア一覧を全条件に適用し、1変異・1,048,576 bp窓・1条件ごとにAPIを1回呼ぶ。[公式の条件選択ガイド](https://www.alphagenomedocs.com/colabs/v2/choosing_a_condition.html)・[複数スコアの取得例](https://www.alphagenomedocs.com/colabs/v2/splicing_variant.html#score-the-variant)

```powershell
.\.venv\Scripts\ag2-variant-scorer.exe score --input examples/apol4.tsv --reference 'C:\path\to\hg38.fa' `
  --condition 'UBERON:0000053_CL:0000573' --condition 'UBERON:0000053_CL:0002586' `
  --scorer SC_GENE_EXPRESSION --scorer RNA_SEQ `
  --scorer SPLICE_SITES --scorer SPLICE_SITE_USAGE --scorer SPLICE_JUNCTIONS `
  --output-dir 'C:\path\to\data\ag2-retina-rpe' --api-key-file 'C:\path\to\private\api_key.txt'
```

`--output-dir` には新しいディレクトリを指定する。`ALPHAGENOME_API_KEY` 環境変数を使う場合は `--api-key-file` を省略できる。この2条件のスプライシングは探索的な補完予測（`imputed_exploratory`）にあたる。解釈時はトラックごとの `imputation_mode` を確認する。[公式の補完予測ガイド](https://www.alphagenomedocs.com/user_guides/imputation.html)

| 出力 | 内容 |
| --- | --- |
| `scores.tsv` | 遺伝子・接合部・トラック別のスコア。`input_id`、`requested_condition`、`requested_scorer`、`raw_score`、`imputation_mode` を保持 |
| `conditions.tsv` | 指定した参照条件のメタデータ |
| `requests.tsv` | 要求ごとの状態とスコア行数 |
| `splicing_summary.tsv` | スプライシング3種をすべて指定した場合に、変異・条件ごとの成分値と統合値を保存 |
| `variants.tsv` | 入力変異の控え |
| `request.json` | 指定条件・スコア、モデルとSDK版などの実行記録 |
| `run.sqlite` | 中断後の再開に使う進捗と取得済み応答 |
| `rejected.tsv` | REF不一致などで除外した入力行と理由。該当時のみ作成 |

入力の不正行は記録して、処理可能な行を続行する。スコア取得中にAPIエラーが起きた場合は、取得済みの結果を保存して処理を打ち切る。どちらの場合も終了コードは `1`。

統合スプライシングスコアは、各成分について有限な `observed` トラックの最大絶対値を求め、`SPLICE_SITES + SPLICE_SITE_USAGE + SPLICE_JUNCTIONS / 5` で計算する。条件ごとに集約し、成分が欠ける場合は統合値を空欄にして `missing_observed_components` に不足成分を記録する。黄斑の錐体とRPEの例では統合値が空欄になり、探索的予測の個別値は `scores.tsv` に残る。[公式の統合スコア](https://www.alphagenomedocs.com/colabs/v2/splicing_variant.html#combine-the-splicing-scores)

## 条件パネルとリクエスト件数の見積もり

`--panel` にTOMLファイルを指定すると、参照条件ごとにscorerの一覧を設定できる。[網膜・RPEの設定例](panels/retina-rpe.toml)は、黄斑と周辺網膜の錐体・RPEの4条件で `SC_GENE_EXPRESSION` を取得する。パネルと `--condition`／`--scorer` による直接指定は実行ごとにどちらかを選ぶ。

```powershell
.\.venv\Scripts\ag2-variant-scorer.exe score --input examples/apol4.tsv --reference 'C:\path\to\hg38.fa' `
  --panel panels/retina-rpe.toml --dry-run
```

`--dry-run` はローカルでREFを照合し、照合を通った入力行数・条件数・予定スコア取得件数を表示する。APIキーなしで実行できる。取得するには `--dry-run` を `--output-dir 'C:\path\to\data\ag2-retina-rpe'` に替え、APIキーを設定する。

スコア取得の予定件数は、REF照合を通った入力行数と条件数の積。1,000行を4条件で評価すると4,000件、20条件では20,000件になる。1条件内の複数scorerは1件にまとめる。実行時の `request.json` には、実際に指定した条件ごとのscorerと使用したパネルを記録する。[AG2の条件指定](https://www.alphagenomedocs.com/user_guides/ag1_vs_ag2.html#conditions)

## 中断からの再開と実行間キャッシュ

`score` は取得のたびにスコアと進捗を出力先の `run.sqlite` に保存する。中断・取得エラーの後は同じ実行ディレクトリを指定して、失敗した取得と未実行の取得を再開する。元の入力と条件は実行ディレクトリから読み込む。取得済みのスコアだけで完了する場合は認証も省略できる。

```powershell
.\.venv\Scripts\ag2-variant-scorer.exe resume --run-dir 'C:\path\to\data\ag2-retina-rpe' `
  --api-key-file 'C:\path\to\private\api_key.txt'
```

別の実行で取得済みスコアを再利用する場合は `--cache-dir 'C:\path\to\cache\ag2'` を追加する。指定は任意で、再開時は保存済みの指定を引き継ぐ。同じ変異・窓・条件・モデル指定・SDK版でscorer設定が一致する結果を再利用し、未取得のscorerを条件ごとにまとめる。行が返らなかった取得は、その実行内の完了として扱う。別実行間のキャッシュにはスコア行がある応答を保存する。

`--dry-run --cache-dir ...` では、既存キャッシュだけを読み、命中したscorer件数と新規取得件数の上限を追加表示する。同じ変異・条件を複数の入力行で指定した場合は、最初の取得結果を同じ実行の後続行でも使えるため、実行件数が上限より少なくなることがある。再開時は、実行DBに保存済みの応答で共有キャッシュの不足分を補う。実行後の `requests.tsv` には `api`／`cache` と、元の取得時刻 `retrieved_at_utc` を記録する。`request.json` の `attempted_score_calls` はCLIが開始した `score_variant()` の累計で、SDK内部の再試行はこの件数に含めない。[SDKの再試行設定](https://github.com/google-deepmind/alphagenome/blob/v2.0.0/src/alphagenome/protos/v2/grpc_service_config.json)

再開には同じモデル指定・SDK・scorer設定を使う。AG2の応答からサービス側のモデル改訂IDは得られないため、更新後の結果を取得したい場合は新しい実行と新しいキャッシュを使う。[AG2の応答仕様](https://github.com/google-deepmind/alphagenome/blob/v2.0.0/src/alphagenome/protos/v2/dna_model_service.proto#L193-L222)

スコア取得中にCtrl+Cで中断すると、取得済みの結果を出力して終了コード130を返す。プロセスが強制終了した場合は、最後にSQLiteへ保存された結果から再開できる。実行DBとキャッシュはローカルのファイルシステムに置く。Linuxの排他と再開はWSL2のext4で確認した。

## 保存済みスコアから遺伝子候補を調べる

`score` の出力ディレクトリを `--run-dir` に渡すと、保存済みのスコアから遺伝子別の集約と変異別の評価範囲を出力する。

```powershell
.\.venv\Scripts\ag2-variant-scorer.exe summarize --run-dir 'C:\path\to\data\ag2-retina-rpe' `
  --output-dir 'C:\path\to\data\ag2-retina-rpe-summary'
```

| 出力 | 内容 |
| --- | --- |
| `gene_summary.tsv` | 変異・遺伝子ID・scorer別の順位、符号付き最小値・最大値とその出典、指定時の `meets_threshold` |
| `variant_summary.tsv` | 変異・scorer別の評価範囲、条件ごとの状態、閾値を指定した場合の候補判定 |
| `summary.json` | 元の実行記録、集約した条件と閾値 |

`SC_GENE_EXPRESSION` と `RNA_SEQ` はそれぞれ集約し、学習時に実測データがある `observed` トラックの有限値を使う。順位は変化の絶対値で付け、閾値なしでも出力する。抽出基準を指定する場合は `--threshold SC_GENE_EXPRESSION=0.2` のようにscorerごとに渡す。`0.2` は指定方法の例。候補なしは、対象条件をすべて評価でき、設定した閾値以上の効果が見つからなかった場合に記録する。取得エラーや未実行などの状態は候補判定と別に残す。

## 入力・出力データの保管

取得時は変異の座標・REF・ALT、予測窓、条件、scorerをAlphaGenome APIへ送る。`input_id` とTSVの追加列はローカルで扱う。`variants.tsv` には追加列も含めた入力の控えが保存される。実行DBには変異と結果が、共有キャッシュにはスコア応答が残る。実データ、認証ファイル、出力先とキャッシュはリポジトリの外に置く。公開例用の `runs/` と `cache/` はGitの除外設定に含まれる。

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

`ag2_variant_scorer/reference.py` のFASTA読み取りとSNVのREF照合は、`atollray/alphagenome-variant-scorer` のコミット `9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc` にある `Reference` と `verified_normalize` の該当部分を基にしている。同コミットのApache-2.0ライセンスを `LICENSE` に収録した。[FASTA読み取りのコピー元](https://github.com/atollray/alphagenome-variant-scorer/blob/9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc/src/alphagenome_variant_scorer/reference.py)・[REF照合のコピー元](https://github.com/atollray/alphagenome-variant-scorer/blob/9b768c199dfcfc4c7e2b3acbb003c8c7c54545dc/src/alphagenome_variant_scorer/variants.py)

## ライセンスとAPI利用条件

このCLIのソースコードにはApache-2.0が適用される。AlphaGenome APIと予測結果にはGoogleの利用条件が適用される。研究用APIと出力の利用範囲、再配布、臨床判断に関する条件は公式案内を参照する。[CLIのライセンス](LICENSE)・[AlphaGenome利用規約](https://alphagenome.google/terms)・[公式FAQ](https://www.alphagenomedocs.com/faqs.html#what-terms-of-use-apply-and-is-commercial-access-available)
