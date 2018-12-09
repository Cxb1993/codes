set terminal latex size 3.5,2.62
set output 'trapezoidal.tex'
set xrange [0:4]
plot "trapezoidal.txt" title 'log(D) VS log(N)' with lines linestyle 1
set output