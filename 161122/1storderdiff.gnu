set terminal epslatex size 13.7cm,13.7cm color colortext

# Line width of the axes
set border linewidth 2

# Axes label
set xlabel '$x$'
set ylabel '$f(x)$'
set format '$%g$'

# Axes ranges
set xrange [*:*]
set yrange [*:*]

# Axes tics
set xtics (0, '$\pi/4$' pi/4, '$\pi/2$' pi/2, '$3\pi/4$' 3*pi/4,\
            '$\pi$' pi)
# set ytics (0, 0.5, 1)
set tics scale 1.5

set key at ((pi) - 0.5),-2

# First order plot for sine function with 10 intervals
set output '1stsin10.tex'
plot '1stsinfwd10.dat' with lines ls 4 title 'Forward',\
'1stsinbwd10.dat' with lines ls 2 title 'Backward',\
'1stsincnt10.dat' with lines ls 3 title 'Central',\
(3*cos(3*x)) with lines ls 1 title 'Analytical'

# First order plot for sine function with 100 intervals
set output '1stsin100.tex'
plot '1stsinfwd100.dat' with lines ls 4 title 'Forward',\
'1stsinbwd100.dat' with lines ls 2 title 'Backward',\
'1stsincnt100.dat' with lines ls 3 title 'Central',\
(3*cos(3*x)) with lines ls 1 title 'Analytical'

set key at ((pi) - 0.2),2.7

# First order plot for cosine function with 10 intervals
set output '1stcos10.tex'
plot '1stcosfwd10.dat' with lines ls 4 title 'Forward',\
'1stcosbwd10.dat' with lines ls 2 title 'Backward',\
'1stcoscnt10.dat' with lines ls 3 title 'Central',\
(-3*sin(3*x)) with lines ls 1 title 'Analytical'

# First order plot for cosine function with 100 intervals
set output '1stcos100.tex'
plot '1stcosfwd100.dat' with lines ls 4 title 'Forward',\
'1stcosbwd100.dat' with lines ls 2 title 'Backward',\
'1stcoscnt100.dat' with lines ls 3 title 'Central',\
(-3*sin(3*x)) with lines ls 1 title 'Analytical'
