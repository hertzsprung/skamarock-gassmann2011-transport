set style data linespoints

set logscale
set format x "10^{%L}"
set format y "10^{%L}"

set xlabel "mean dx"

set xrange [1e-1:5e-4]
set yrange [1e-10:1]

set key outside top center

set multiplot layout 1,2

set ylabel "l2 error"

plot 'results/sk-uniform.dat' using 1:2 lc 1 dt 1 title "skam-gass uniform", \
     'results/sk-nonuniform.dat' using 1:2 lc 1 dt 3 title "skam-gass nonuniform", \
     'results/cf-uniform.dat' using 1:2 lc 2 dt 1 title "cubicFit uniform", \
     'results/cf-nonuniform.dat' using 1:2 lc 2 dt 3 title "cubicFit nonuniform", \
     'results/cf-corr-uniform.dat' using 1:2 lc 3 dt 1 title "cubicFit corr uniform", \
     'results/cf-corr-nonuniform.dat' using 1:2 lc 3 dt 3 title "cubicFit corr nonuniform", \
     x * 5e-1 lc rgbcolor "black" dt 1 lw 2 title '1st order', \
     x**2 * 4 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**3 * 1e2 lc rgbcolor "black" dt 3 lw 2 title '3rd order', \
     x**4 * 1e3 lc rgbcolor "black" dt 4 lw 2 title '4th order'

set ylabel "linfty error"

plot 'results/sk-uniform.dat' using 1:3 lc 1 dt 1 title "skam-gass uniform", \
     'results/sk-nonuniform.dat' using 1:3 lc 1 dt 3 title "skam-gass nonuniform", \
     'results/cf-uniform.dat' using 1:3 lc 2 dt 1 title "cubicFit uniform", \
     'results/cf-nonuniform.dat' using 1:3 lc 2 dt 3 title "cubicFit nonuniform", \
     'results/cf-corr-uniform.dat' using 1:3 lc 3 dt 1 title "cubicFit corr uniform", \
     'results/cf-corr-nonuniform.dat' using 1:3 lc 3 dt 3 title "cubicFit corr nonuniform", \
     x * 5e-1 lc rgbcolor "black" dt 1 lw 2 title '1st order', \
     x**2 * 4 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**3 * 1e2 lc rgbcolor "black" dt 3 lw 2 title '3rd order', \
     x**4 * 1e3 lc rgbcolor "black" dt 4 lw 2 title '4th order'
