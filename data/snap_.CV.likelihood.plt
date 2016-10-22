#
# hold-out likelihood (Sat Oct 22 13:31:11 2016)
#

set title "hold-out likelihood"
set key bottom right
set autoscale
set grid
set xlabel "communities"
set ylabel "likelihood"
set tics scale 2
set terminal png size 1000,800
set output 'snap_.CV.likelihood.png'
plot 	"snap_.CV.likelihood.tab" using 1:2 title "" with linespoints pt 6
