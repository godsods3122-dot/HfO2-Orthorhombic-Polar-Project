set terminal pdf enhanced color font ",24"
set palette defined ( 0  "green", 5 "yellow", 10 "red" )
set output 'bulkek.pdf' 
set style data linespoints
unset key
set pointsize 0.8
#set xtics font ",24"
#set ytics font ",24"
#set ylabel font ",24"
set ylabel offset 0.5,0
set xrange [0:    5.67094]
emin=   -0.582571
emax=   25.318063
set yrange [emin: emax]
set ylabel "Frequency (THz)"
set xtics ("G  "    0.00000,"X  "    0.61927,"S  "    1.26409,"Y  "    1.88336,"G  "    2.52818,"Z  "    3.14276,"U  "    3.76203,"R  "    4.40685,"T  "    5.02612,"Z  "    5.67094)
set arrow from    0.61927, emin to    0.61927, emax nohead
set arrow from    1.26409, emin to    1.26409, emax nohead
set arrow from    1.88336, emin to    1.88336, emax nohead
set arrow from    2.52818, emin to    2.52818, emax nohead
set arrow from    3.14276, emin to    3.14276, emax nohead
set arrow from    3.76203, emin to    3.76203, emax nohead
set arrow from    4.40685, emin to    4.40685, emax nohead
set arrow from    5.02612, emin to    5.02612, emax nohead
# please comment the following lines to plot the fatband 
plot 'bulkek.dat' u 1:2  w lp lw 0.1 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2
 
# uncomment the following lines to plot the fatband 
#plot 'bulkek.dat' u 1:2:3  w lp lw 2 pt 7  ps 0.2 lc palette, 0 w l lw 2
# uncomment the following lines to plot the spin if necessary
#plot 'bulkek.dat' u 1:2 w lp lw 2 pt 7  ps 0.2, \
     'bulkek.dat' u 1:2:($3/6):($4/6) w vec
