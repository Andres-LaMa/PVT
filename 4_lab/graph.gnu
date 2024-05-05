#!/usr/bin/gnuplot
set termoption enhanced
set terminal postscript eps enhanced color solid font "Arial, 16"
# set terminal png size 800, 600
set output "speedup.ps"
set style line 1 lc rgb "0xDC143C" lt 1 lw 4 pt 9 ps 1
set style line 2 lc rgb "green" lt 1 lw 4 pt 9 ps 1
set style line 3 lc rgb "purple" lt 1 lw 4 pt 9 ps 1
set style line 4 lc rgb "red" lt 1 lw 4 pt 9 ps 1
set style line 5 lc rgb "dark-blue" lt 1 lw 4 pt 9 ps 1
set style line 6 lc rgb "black" lt 1 lw 4 pt 9 ps 1

set border lw 2
set grid

set key top left
set xlabel "Количество потоков"
set ylabel "Коэффициент ускорения" rotate by 90

set xtics 1
set mxtics

set xrange [0:8]
set yrange [0:2]

set format x "%6.0f"
set format y "%.6f"

plot x title "line speedup^a" with linespoints ls 1, \
"cri.csv" using 1:2 title "critical" with linespoints ls 2, \
"ato.csv" using 1:2 title "atomic" with linespoints ls 3, \
"n_block.csv" using 1:2 title "n blocked" with linespoints ls 4, \
"extra_calc.csv" using 1:2 title "extra calc" with linespoints ls 5, \
"extra_mem.csv" using 1:2 title "extra mem" with linespoints ls 6