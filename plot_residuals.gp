set terminal x11 size 1024,768

set grid
set title 'OpenCFD-EC Residuals History'
set xlabel 'Iteration Step (Kstep)'
set ylabel 'RMS Residual (Log Scale)'

set logscale y
set format y "10^{%L}"

set style data lines

plot 'Residual.dat' using 1:7 title 'RMS \rho' linewidth 1.5, \
     'Residual.dat' using 1:8 title 'RMS \rho u' linewidth 1.5, \
     'Residual.dat' using 1:9 title 'RMS \rho v' linewidth 1.5, \
     'Residual.dat' using 1:10 title 'RMS \rho w' linewidth 1.5, \
     'Residual.dat' using 1:11 title 'RMS \rho E' linewidth 1.5

while (1) {
    pause 2       # Wait for 2 seconds
    replot        # Refresh the screen with new data
}
