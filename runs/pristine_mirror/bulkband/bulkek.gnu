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
set xrange [0:    5.55602]
emin=   -0.572491
emax=   24.111575
set yrange [emin: emax]
set ylabel "Frequency (THz)"
set xtics ("G  "    0.00000,"X  "    0.60380,"S  "    1.23251,"Y  "    1.83632,"G  "    2.46503,"Z  "    3.09099,"U  "    3.69480,"R  "    4.32351,"T  "    4.92731,"Z  "    5.55602)
set arrow from    0.60380, emin to    0.60380, emax nohead
set arrow from    1.23251, emin to    1.23251, emax nohead
set arrow from    1.83632, emin to    1.83632, emax nohead
set arrow from    2.46503, emin to    2.46503, emax nohead
set arrow from    3.09099, emin to    3.09099, emax nohead
set arrow from    3.69480, emin to    3.69480, emax nohead
set arrow from    4.32351, emin to    4.32351, emax nohead
set arrow from    4.92731, emin to    4.92731, emax nohead
# please comment the following lines to plot the fatband 
plot 'bulkek.dat' u 1:2  w lp lw 0.1 pt 7  ps 0.2 lc rgb 'black', 0 w l lw 2
 
# uncomment the following lines to plot the fatband 
#plot 'bulkek.dat' u 1:2:3  w lp lw 2 pt 7  ps 0.2 lc palette, 0 w l lw 2
# uncomment the following lines to plot the spin if necessary
#plot 'bulkek.dat' u 1:2 w lp lw 2 pt 7  ps 0.2, \
     'bulkek.dat' u 1:2:($3/6):($4/6) w vec
