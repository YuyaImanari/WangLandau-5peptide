# WL_5peptide

Wang Landau アルゴリズムを用いて、5 つのアミノ酸からなるペプチドのエネルギー分布を予測するプログラムです。

## 環境

- pipenv : version 2022.5.2
- python : version 3.9.6

## 準備

1. リポジトリをクローンする

1. パッケージをインストールする

    ```bash
    # pipenv を使用している場合
    pipenv install
    # pipenv install -r ./requirements.txt でも可

    # pipenv を使用していない場合
    pip install -r requirements.txt
    ```

## 実行

コマンドライン引数 `N` が ステップ数 (numiter) に相当します。

```bash
# pipenv を使用している場合
pipenv run python main.py N

# シェルに入ってから、★のように実行しても可
pipenv shell

# pipenv を使用していない場合 ★
python main.py N
```
