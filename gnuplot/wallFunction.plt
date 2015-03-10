#!/bin/gnuplot
set terminal pdf
set output 'wallFunction.pdf'
E = 9.1
k = 0.42
ak= 0.01
bk= 5
set grid
set contour
set cntrparam levels discrete 0
unset surface
set view map
set isosamples 1000,1000

set key bottom
set box
set xrange[0.1:500]
set yrange[0.0:30]
set log x
#set log y

set xlabel "y+"
set ylabel "u+"

#Spalding Wall function
fs(x,y) = y - x + (1.0/E)*(exp(k*y) - 1 - k*y- 0.5*(k*y)**2 - (1.0/6.0)*(k*y)**3)

#Enhanced Wall function
fk(x,y) = exp((-ak*(x)**4)/(1.0+bk*x)) * x + exp(1.0/((-ak*(x)**4)/(1.0+bk*x)))*(1.0/k)*log(E*x) - y                                   

#Standard Wall function
fw(x,y) = x <= 11.81 ? x - y : (1.0/k)*log(E*x) - y

set title "Wall function"


splot  fs(x,y) ls 1 title "Spalding Wall function"\
      ,fk(x,y) ls 2 title "Enhanced Wall function"\
      ,fw(x,y) ls 3 title "Standard Wall function"


pause -1
