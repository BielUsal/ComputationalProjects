set view map
set xlabel "x (m)"
set ylabel "y (m)"
splot "V-SOR.dat" u 1:2:3 ps 4 pt 5 palette
pause mouse
