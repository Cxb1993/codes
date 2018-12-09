set terminal epslatex size 13.7cm,13.7cm color colortext
set output 'xygrid.tex'

# Line width of the axes
set border linewidth 2

# Axes label
set xlabel '$x$'
set ylabel '$y$'
set format '$%g$'

# Axes ranges
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]

# Axes tics
set xtics (0, 0.5, 1)
set ytics (0, 0.5, 1)
set tics scale 1.5

# Plot
plot 'xygrid.dat' with lines notitle