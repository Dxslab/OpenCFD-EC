# Gnuplot script to plot residuals from Residual.dat
# Run this script using: gnuplot plot_residuals.gp

# Set terminal to screen for real-time display
set terminal x11 size 1024,768 enhanced font 'Verdana,10' persist
# set output 'Residuals.png' # Comment out output to file

# Set grid and labels
set grid
set title 'OpenCFD-EC Residuals History' font ',14'
set xlabel 'Iteration Step (Kstep)' font ',12'
set ylabel 'RMS Residual (Log Scale)' font ',12'

# Use logarithmic scale for Y axis
set logscale y
set format y "10^{%L}"

# Configure plot style
set style data lines
set key outside right top

# In Residual.dat:
# Column 1 is Kstep
# Columns 2 to (1+NVAR) are Res_max
# Columns (2+NVAR) to end are Res_rms
# Usually NVAR is 5 (Navier-Stokes) or 6/7 (with turbulence models)

# Assuming NVAR=5, RMS columns start at 7
# (Kstep, Max(1..5), RMS(1..5)) => Kstep is col 1, RMS1 is col 7, RMS2=8, RMS3=9, RMS4=10, RMS5=11
# We will plot RMS residuals. If your NVAR depends on the turbulence model, adjust the column numbers below.

plot 'Residual.dat' using 1:7 title 'RMS \rho' linewidth 1.5, \
     'Residual.dat' using 1:8 title 'RMS \rho u' linewidth 1.5, \
     'Residual.dat' using 1:9 title 'RMS \rho v' linewidth 1.5, \
     'Residual.dat' using 1:10 title 'RMS \rho w' linewidth 1.5, \
     'Residual.dat' using 1:11 title 'RMS \rho E' linewidth 1.5

# Real-time update loop
while (1) {
    pause 2       # Wait for 2 seconds
    replot        # Refresh the screen with new data
}
