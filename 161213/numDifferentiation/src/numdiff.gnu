set terminal epslatex size 13.7cm,13.7cm color colortext

set datafile separator "\t"

# Line width of the axes
set border linewidth 2

unset logscale

# Axes label
set xlabel '$x$'
set ylabel '$f(x)$'
set format '$%g$'

# Axes ranges

set xrange [0:pi]
set yrange [*:*]

# Axes tics
set xtics ('$0$' 0, '$\pi/4$' pi/4, '$\pi/2$' pi/2, '$3\pi/4$' 3*pi/4, '$\pi$' pi)
set tics scale 1.5

set key at ((pi) - 0.5),-2

# First derivative plot for a sine function with 10 intervals
set output 'D:/PhD/Tasks/004/report/firstsin10.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/firstsin10.dat' using 1:2 with lines ls 4 title "Analytical",\
'D:/PhD/Tasks/004/numDifferentiation/data/firstsin10.dat' using 1:3 with lines ls 1 title "Forward",\
'D:/PhD/Tasks/004/numDifferentiation/data/firstsin10.dat' using 1:5 with lines ls 2 title "Backward",\
'D:/PhD/Tasks/004/numDifferentiation/data/firstsin10.dat' using 1:7 with lines ls 3 title "Central"

# First derivative plot for a sine function with 100 intervals
set output 'D:/PhD/Tasks/004/report/firstsin100.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/firstsin100.dat' using 1:2 with lines ls 4 title "Analytical",\
'D:/PhD/Tasks/004/numDifferentiation/data/firstsin100.dat' using 1:3 with lines ls 1 title "Forward",\
'D:/PhD/Tasks/004/numDifferentiation/data/firstsin100.dat' using 1:5 with lines ls 2 title "Backward",\
'D:/PhD/Tasks/004/numDifferentiation/data/firstsin100.dat' using 1:7 with lines ls 3 title "Central"

set key at ((pi) - 0.1),9

# Second derivative plot for a sine function with 10 intervals
set output 'D:/PhD/Tasks/004/report/secondsin10.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/secondsin10.dat' using 1:2 with lines ls 4 title "Analytical",\
'D:/PhD/Tasks/004/numDifferentiation/data/secondsin10.dat' using 1:3 with lines ls 1 title "Forward",\
'D:/PhD/Tasks/004/numDifferentiation/data/secondsin10.dat' using 1:5 with lines ls 2 title "Backward",\
'D:/PhD/Tasks/004/numDifferentiation/data/secondsin10.dat' using 1:7 with lines ls 3 title "Central"

# Second derivative plot for a sine function with 100 intervals
set output 'D:/PhD/Tasks/004/report/secondsin100.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/secondsin100.dat' using 1:2 with lines ls 4 title "Analytical",\
'D:/PhD/Tasks/004/numDifferentiation/data/secondsin100.dat' using 1:3 with lines ls 1 title "Forward",\
'D:/PhD/Tasks/004/numDifferentiation/data/secondsin100.dat' using 1:5 with lines ls 2 title "Backward",\
'D:/PhD/Tasks/004/numDifferentiation/data/secondsin100.dat' using 1:7 with lines ls 3 title "Central"

# Axes label
set xlabel 'log$_{10}N$'
set ylabel 'log$_{10}E$'
set format '$%g$'

# Axes ranges

unset xrange
unset yrange
set xrange [*:*]
set yrange [*:*]

# Axes tics
unset xtics
set tics scale 1.5

unset key

set logscale x;

# Error of first derivative plot for a sine function for forward difference method
set output 'D:/PhD/Tasks/004/report/firstsinfwderror.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/firstsinerror.dat' using 1:5 with lines ls 1

# Error of first derivative plot for a sine function for backward difference method
set output 'D:/PhD/Tasks/004/report/firstsinbwderror.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/firstsinerror.dat' using 1:8 with lines ls 1

# Error of first derivative plot for a sine function for central difference method
set output 'D:/PhD/Tasks/004/report/firstsincnterror.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/firstsinerror.dat' using 1:11 with lines ls 1

# Error of second derivative plot for a sine function for forward difference method
set output 'D:/PhD/Tasks/004/report/secondsinfwderror.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/secondsinerror.dat' using 1:5 with lines ls 1

# Error of second derivative plot for a sine function for backward difference method
set output 'D:/PhD/Tasks/004/report/secondsinbwderror.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/secondsinerror.dat' using 1:8 with lines ls 1

# Error of second derivative plot for a sine function for central difference method
set output 'D:/PhD/Tasks/004/report/secondsincnterror.tex'
plot 'D:/PhD/Tasks/004/numDifferentiation/data/secondsinerror.dat' using 1:11 with lines ls 1
