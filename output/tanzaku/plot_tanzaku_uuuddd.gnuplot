# === 共通設定 ===
set datafile separator ","
set title "Energy vs. Berry Curvature"
set xlabel "Energy"
set ylabel "Values"
set grid
set key outside
set xrange [-0.3:0.3]
set yrange [-1:1]
set terminal pngcairo size 800,600

# === ファイル名・ラベル・色の定義 ===
files = "./data/uuuddd_lambda0p1_j0p25_mu0_100_300_0.2.dat \
         ./data/uuuddd_lambda0p1_j0p25_mu0_200_300_0.2.dat \
         ./data/uuuddd_lambda0p1_j0p25_mu0_300_300_0.2.dat \
         ./data/uuuddd_lambda0p1_j0p25_mu0_400_300_0.2.dat \
         ./data/uuuddd_lambda0p1_j0p25_mu0_500_300_0.2.dat \
         ./data/uuuddd_lambda0p1_j0p25_mu0_600_300_0.2.dat"

labels = "100 200 300 400 500 600"
colors = "#8c564b #d62728 #1f77b4 #2ca02c #ff7f0e #9467bd"  # 茶, 赤, 青, 緑, オレンジ, 紫

# === bcd_y のプロット ===
set output "plot_uuuddd_only_bcd_y_compare.png"
plot for [i=1:6] word(files, i) using 1:4 with lines \
     linecolor rgb sprintf("%s", word(colors, i)) \
     title sprintf("bcd y %s", word(labels, i))

# === bcd_x のプロット ===
set output "plot_uuuddd_only_bcd_x_compare.png"
plot for [i=1:6] word(files, i) using 1:3 with lines \
     linecolor rgb sprintf("%s", word(colors, i)) \
     title sprintf("bcd x %s", word(labels, i))

# === bc のプロット ===
set output "plot_uuuddd_only_bc_compare.png"
plot for [i=1:6] word(files, i) using 1:2 with lines \
     linecolor rgb sprintf("%s", word(colors, i)) \
     title sprintf("bc %s", word(labels, i))
