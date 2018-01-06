set pm3d map
set terminal png
set view map
set title "Local density in a nuclear system"

unset cbtics
set palette rgbformulae -5, 3, -7
set cblabel "rho-value"
set xlabel "R"
set ylabel "Z"
DEBUG_TERM_HTIC = 75
DEBUG_TERM_VTIC = 75
set output "result.png"
splot "result.dat" matrix
