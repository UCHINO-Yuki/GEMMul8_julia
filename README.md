# GEMMul8_Julia

## 1. Julia と作業フォルダを確認する

### Windows PowerShell

```powershell
Set-Location "C:\path\to\ozaki-julia"
julia --version
```

### Linux / macOS

```bash
cd "/path/to/ozaki-julia"
julia --version
```

フォルダのパスは実際の保存先に置き換えてください。
`julia` が見つからない場合は Julia をインストールし、`julia` コマンドのあるディレクトリを PATH に追加してから、ターミナルを開き直してください。

## 2. コードを実行する

PowerShell・Linux・macOS・WSL で共通です。

```text
julia --project=. test.jl
```

実行が完了すると、`test.jl` と同じフォルダの `results.csv` を生成・上書きします。
現在の設定は `m = 8`、`n = 8`、`k = 8192`、`phi = 0, 1, 2, 4, 8, 12`、法の個数 `3:40` です。
変更する場合は `test.jl` 内の対応する値を編集してください。

`test.jl` は `MultiFloats`・`CSV`・`DataFrames` を確認し、見つからないパッケージだけを自動追加します。
`gemmul8.jl` が直接読み込む `LinearAlgebra` は Julia の標準ライブラリです。

## 3. 図を作成する

```text
julia --project=. plot.jl
```

`results.csv` を描画する。
`plot.jl` は `CSV`・`DataFrames`・`Plots`・`LaTeXStrings` を確認し、見つからないものだけを自動追加します。

出力先は `plot.jl` と同じフォルダです。

| 出力 | 内容 |
| --- | --- |
| `plot_S.png` | Float32 入力：FP64 CRT・Exact CRT・SGEMM の比較 |
| `plot_D.png` | Float64 入力：Pseudo-DD CRT・Exact CRT・DGEMM の比較 |

### パッケージの追加先

`--project=.` は、カレントディレクトリを有効な Julia プロジェクトとして指定します。
自動追加が必要な場合、`Pkg.add` はそのプロジェクトの `Project.toml` と `Manifest.toml` を作成・更新します。
各スクリプトは `Base.find_package` で検索するため、Julia の読み込みパス上ですでに見つかるパッケージは追加しません。
初回の追加・コンパイルには時間がかかることがあります。
