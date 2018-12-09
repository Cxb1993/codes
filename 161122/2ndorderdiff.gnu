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

set key default
set key at ((pi) - 0.1),8.7

# First order plot for sine function with 10 intervals
set output '2ndsin10.tex'
plot '2ndsinfwd10.dat' with lines ls 4 title 'Forward',\
'2ndsinbwd10.dat' with lines ls 2 title 'Backward',\
'2ndsincnt10.dat' with lines ls 3 title 'Central',\
(-9*sin(3*x)) with lines ls 1 title 'Analytical'

# First order plot for sine function with 100 intervals
set output '2ndsin100.tex'
plot '2ndsinfwd100.dat' with lines ls 4 title 'Forward',\
'2ndsinbwd100.dat' with lines ls 2 title 'Backward',\
'2ndsincnt100.dat' with lines ls 3 title 'Central',\
(-9*sin(3*x)) with lines ls 1 title 'Analytical'

set key at ((pi) - 0.4),8.7

# First order plot for cosine function with 10 intervals
set output '2ndcos10.tex'
plot '2ndcosfwd10.dat' with lines ls 4 title 'Forward',\
'2ndcosbwd10.dat' with lines ls 2 title 'Backward',\
'2ndcoscnt10.dat' with lines ls 3 title 'Central',\
(-9*cos(3*x)) with lines ls 1 title 'Analytical'

# First order plot for cosine function with 100 intervals
set output '2ndcos100.tex'
plot '2ndcosfwd100.dat' with lines ls 4 title 'Forward',\
'2ndcosbwd100.dat' with lines ls 2 title 'Backward',\
'2ndcoscnt100.dat' with lines ls 3 title 'Central',\
(-9*cos(3*x)) with lines ls 1 title 'Analytical'
