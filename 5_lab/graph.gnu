#!/usr/bin/gnuplot
set termoption enhanced
set terminal postscript eps enhanced color solid font "Arial, 16"
# set terminal png size 800, 600

set style line 1 lc rgb "dark-blue" lt 1 lw 4 pt 9 ps 1
set style line 2 lc rgb "green" lt 1 lw 4 pt 9 ps 1
set style line 3 lc rgb "purple" lt 1 lw 4 pt 9 ps 1
set style line 4 lc rgb "red" lt 1 lw 4 pt 9 ps 1

set border lw 2
set grid

set key top left
set xlabel "Количество потоков"
set ylabel "Коэффициент ускорения" rotate by 90

set xtics 1
set mxtics

set xrange [0:8]
set yrange [0:4]

set format x "%6.0f"
set format y "%.1f"


set output "speedup_1000.ps"
plot x title "line speedup" with linespoints ls 1, \
"quik_sort_th=1000.csv" using (strlen(strcol(1)) == 8 ? $2 : NaN):3 title "1000000" with linespoints ls 2, \
"quik_sort_th=1000.csv" using (strlen(strcol(1)) == 9 ? $2 : NaN):3 title "10000000" with linespoints ls 3, \
"quik_sort_th=1000.csv" using (strlen(strcol(1)) == 10 ? $2 : NaN):3  title "100000000" with linespoints ls 4

set output "speedup_10000.ps"
plot x title "line speedup" with linespoints ls 1, \
"quik_sort_th=10000.csv" using (strlen(strcol(1)) == 8 ? $2 : NaN):3 title "1000000" with linespoints ls 2, \
"quik_sort_th=10000.csv" using (strlen(strcol(1)) == 9 ? $2 : NaN):3 title "10000000" with linespoints ls 3, \
"quik_sort_th=10000.csv" using (strlen(strcol(1)) == 10 ? $2 : NaN):3  title "100000000" with linespoints ls 4

set output "speedup_100000.ps"
plot x title "line speedup" with linespoints ls 1, \
"quik_sort_th=100000.csv" using (strlen(strcol(1)) == 8 ? $2 : NaN):3 title "1000000" with linespoints ls 2, \
"quik_sort_th=100000.csv" using (strlen(strcol(1)) == 9 ? $2 : NaN):3 title "10000000" with linespoints ls 3, \
"quik_sort_th=100000.csv" using (strlen(strcol(1)) == 10 ? $2 : NaN):3  title "100000000" with linespoints ls 4
