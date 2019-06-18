#-------------------------------------------------------------------------------
# gnuplotの設定
#-------------------------------------------------------------------------------
reset
set nokey                # 凡例の非表示
set xrange [0:1]         # x軸方向の範囲の設定
set yrange [0:1]         # y軸方向の範囲の設定
set size square          # 図を正方形にする

set term gif animate     # 出力をgifアニメに設定
set output "output.gif"  # 出力ファイル名の設定

#-------------------------------------------------------------------------------
# 変数の設定
#-------------------------------------------------------------------------------
n0 = 0    # ループ変数の初期値
n1 = 250  # ループ変数の最大値
dn = 5    # ループ変数の増加間隔

#-------------------------------------------------------------------------------
# ループの開始
#-------------------------------------------------------------------------------
load "loop.plt"