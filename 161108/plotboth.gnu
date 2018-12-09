set output 'both.png
set xrange [0:4]
set xlabel "log(N)"
set ylabel "log(D)"
set title "log(D) VS log(N)"
plot "trapezoidal.txt" title 'Trapezoidal rule' with lines, "simpsons13.txt" title 'Simpsons 1/3 rule' with lines linestyle 2