set terminal epslatex size 13.7cm,13.7cm color colortext

set datafile separator "\t"

# Line width of the axes
set border linewidth 2

# Axes label
set xlabel '$x$'
set ylabel 'Temperature ($^{\circ}$C)'
set format '$%g$'

# Axes ranges

set xrange [*:*]
set yrange [*:*]

# Axes tics
#set xtics ('$0$' 0, '$\pi/4$' pi/4, '$\pi/2$' pi/2, '$3\pi/4$' 3*pi/4, '$\pi$' pi)
set tics scale 1.5

set key at 0.5, 49

# First derivative plot for a sine function with 10 intervals
set output 'D:/PhD/Tasks/005/report161215/ht1dn190t100.tex'
plot 'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t100.dat' using 1:2 with lines ls 1 lw 4 title "t = 0.1",\
'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t1000.dat' using 1:2 with lines ls 2 lw 4 title "t = 1",\
'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t10000.dat' using 1:2 with lines ls 3 lw 4 title "t = 10",\
'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t100000.dat' using 1:2 with lines ls 4 lw 4 title "t = 100",\
'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t1000000.dat' using 1:2 with lines ls 1 lw 4 title "t = 1000",\
'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t2100000.dat' using 1:2 with lines ls 2 lw 4 title "t = 2100",\
'D:\PhD\Tasks\005\heatTransfer1D\data\ht1dn190t3000000.dat' using 1:2 with lines ls 3 lw 4 title "t = 3000"
