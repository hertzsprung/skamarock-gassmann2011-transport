set terminal svg size 800,600 fname 'Verdana' fsize 14
set output 'build/convergence.svg'
set style data linespoints

set logscale
set format x "10^{%L}"
set format y "10^{%L}"

set xlabel "mean Δx"

set xrange [1e-1:5e-4]
set yrange [1e-10:1]

set ylabel "𝓁_2 error"

plot 'build/cf-uniform.dat' using 1:2 lw 2 lc 2 dt 1 title "Uncorrected, uniform mesh", \
     'build/cf-nonuniform.dat' using 1:2 lw 2lc 2 dt 3 title "Uncorrected, nonuniform mesh", \
     'build/cf-corr-uniform.dat' using 1:2 lw 2 lc 3 dt 1 title "Corrected, uniform mesh", \
     'build/cf-corr-nonuniform.dat' using 1:2 lw 2 lc 3 dt 3 title "Corrected, nonuniform mesh", \
     x**2 * 4 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**4 * 1e3 lc rgbcolor "black" dt 4 lw 2 title '4th order'
