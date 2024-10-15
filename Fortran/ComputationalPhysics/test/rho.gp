set view map
set xlabel "x (m)"
set ylabel "y (m)"
splot "rho2.1.dat" u 1:2:3 ps 4 pt 5 palette
pause mouse
