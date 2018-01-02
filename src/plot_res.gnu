set term png
set output 'plot_res.png'
unset key
set view map 

unset cbtics

set xrange [-1:100]
set yrange [-1:100]
set title "test"
set cblabel "Score" 
set cbrange [ 0.14000 : 0.150 ]
set palette rgbformulae -5, 3, -7
#DEBUG_TERM_HTIC = 75
#DEBUG_TERM_VTIC = 75
splot "result.dat" matrix with image