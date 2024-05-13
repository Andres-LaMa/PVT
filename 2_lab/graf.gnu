set termoption enhanced
set terminal postscript eps enhanced color solid font "Arial, 16"
# set terminal png size 800, 600

set style line 1 lc rgb "dark-blue" lt 1 lw 4 pt 9 ps 1
set style line 2 lc rgb "green" lt 1 lw 4 pt 9 ps 1
set style line 3 lc rgb "purple" lt 1 lw 4 pt 9 ps 1

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
set format y "%.1f"

set output "graf.ps"
plot x title "line speedup" with linespoints ls 1, \
"data" using ($1 == 15000 ? $2 : NaN):3 title "15000" with linespoints ls 2, \
"data" using ($1 == 20000 ? $2 : NaN):3 title "20000" with linespoints ls 3, \
"data" using ($1 == 25000 ? $2 : NaN):3 title "25000" with linespoints ls 4