set encoding iso_8859_1
set terminal png truecolor enhanced font ", 60" size 1920, 1680
set output 'wanniercenter3D_Weyl_1.png'
unset key 
set border lw 3 
set xtics offset 0, 0.2
set xtics ("South" 0, "North" 1)
set xtics format "%4.1f" nomirror in
set xlabel "k" 
set xlabel offset 0, 0.7 
set ytics 1.0 
set ytics format "%4.1f" nomirror in
set ylabel "WCC"
set title "Weyl ( 0.44142, 0.50000, 0.01226) "
set ylabel offset 2, 0.0 
set xrange [0: 1.0]
set yrange [0:1]
plot 'wanniercenter3D_Weyl.dat' u 1:    2 w p ps 1.5 pt 7 lw 10 lc rgb '#696969'
