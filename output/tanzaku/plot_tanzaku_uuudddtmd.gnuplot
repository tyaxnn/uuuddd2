# ファイル名と共通設定
set datafile separator ","
set title "Energy vs. Berry Curvature"
set xlabel "Energy"
set ylabel "Values"
set grid
set key outside
set xrange [-3:3]
set yrange [-3:3]

# bcd_y のプロット
set terminal pngcairo size 800,600
set output "plot_uuudddtmd_only_bcd_y_compare.png"
plot "uuudddtmd_141_300_0p6.dat" using 1:4 with lines linecolor rgb "#1f77b4" title "bcd y 141",\
    "uuudddtmd_300_300_0p6.dat" using 1:4 with lines linecolor rgb "#d62728" title "bcd y 300",\
    "uuudddtmd_640_300_0p6.dat" using 1:4 with lines linecolor rgb "#2ca02c" title "bcd y 640"

# bcd_x のプロット
set terminal pngcairo size 800,600
set output "plot_uuudddtmd_only_bcd_x_compare.png"
plot "uuudddtmd_141_300_0p6.dat" using 1:3 with lines linecolor rgb "#2ca02c" title "bcd x 141",\
    "uuudddtmd_300_300_0p6.dat" using 1:3 with lines linecolor rgb "#1f77b4" title "bcd x 300",\
    "uuudddtmd_640_300_0p6.dat" using 1:3 with lines linecolor rgb "#d62728" title "bcd x 640"

# bc のプロット
set terminal pngcairo size 800,600
set output "plot_uuudddtmd_only_bc_compare.png"
plot "uuudddtmd_141_300_0p6.dat" using 1:2 with lines linecolor rgb "#d62728" title "bc 141",\
    "uuudddtmd_300_300_0p6.dat" using 1:2 with lines linecolor rgb "#2ca02c" title "bc 300",\
    "uuudddtmd_640_300_0p6.dat" using 1:2 with lines linecolor rgb "#1f77b4" title "bc 640"
