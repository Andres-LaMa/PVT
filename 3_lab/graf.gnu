#!/usr/bin/gnuplot
set termoption enhanced
set terminal postscript eps enhanced color solid font "Arial, 16"
set output "runge.ps"
set style line 1 lc rgb "0xDC143C" lt 1 lw 4 pt 9 ps 1
set style line 2 lc rgb "green" lt 1 lw 4 pt 9 ps 1
set border lw 2
set grid
set key top left
set xlabel "Количество потоков"
set ylabel "Коэффициент ускорения" rotate by 90
set xtics 1
set mxtics
set xrange [0:8]
set yrange [0:8]
set format x "%6.0f"
set format y "%.6f"
plot "graf.dat" using 1:2 title "runge" with linespoints ls 1, \
"graf.dat" using 3:4 title "line" with linespoints ls 2