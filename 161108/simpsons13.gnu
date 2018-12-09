set terminal latex size 3.5,2.62
set output 'simpsons13.tex'
set xrange [0:4]
plot "simpsons13.txt" title 'log(D) VS log(N)' with lines linestyle 1
set output