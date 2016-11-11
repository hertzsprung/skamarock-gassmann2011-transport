set style data linespoints

set logscale
set format x "10^{%L}"
set format y "10^{%L}"

set xlabel "mean dx"

set xrange [1e-1:1e-3]
set yrange [1e-6:1]

set key outside top center

set multiplot layout 1,2

set ylabel "l2 error"

plot 'results/sk-uniform.dat' using 1:2 lc 1 dt 1 title "skam-gass uniform", \
     'results/sk-nonuniform.dat' using 1:2 lc 1 dt 3 title "skam-gass nonuniform", \
     'results/ls-uniform.dat' using 1:2 lc 2 dt 1 title "ls uniform", \
     'results/ls-nonuniform.dat' using 1:2 lc 2 dt 3 title "ls nonuniform", \
     'results/cf-uniform.dat' using 1:2 lc 3 dt 1 title "cubicFit uniform", \
     'results/cf-nonuniform.dat' using 1:2 lc 3 dt 3 title "cubicFit nonuniform", \
     'results/centred-uniform.dat' using 1:2 lc 4 dt 1 title "centred uniform", \
     'results/centred-nonuniform.dat' using 1:2 lc 4 dt 3 title "centred nonuniform", \
     'results/sknuc-uniform.dat' using 1:2 lc 5 dt 1 title "skam-gass-nonu-centring uniform", \
     'results/sknuc-nonuniform.dat' using 1:2 lc 5 dt 3 title "skam-gass-nonu-centring nonuniform", \
     x * 5e-1 lc rgbcolor "black" dt 1 lw 2 title '1st order', \
     x**2 * 1e1 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**3 * 1e2 lc rgbcolor "black" dt 3 lw 2 title '3rd order'

set ylabel "linfty error"

plot 'results/sk-uniform.dat' using 1:3 lc 1 dt 2 title "skam-gass uniform", \
     'results/sk-nonuniform.dat' using 1:3 lc 1 dt 4 title "skam-gass nonuniform", \
     'results/ls-uniform.dat' using 1:3 lc 2 dt 2 title "ls uniform", \
     'results/ls-nonuniform.dat' using 1:3 lc 2 dt 4 title "ls nonuniform", \
     'results/cf-uniform.dat' using 1:3 lc 3 dt 2 title "cubicFit uniform", \
     'results/cf-nonuniform.dat' using 1:3 lc 3 dt 4 title "cubicFit nonuniform", \
     'results/centred-uniform.dat' using 1:3 lc 4 dt 2 title "centred uniform", \
     'results/centred-nonuniform.dat' using 1:3 lc 4 dt 4 title "centred nonuniform", \
     'results/sknuc-uniform.dat' using 1:3 lc 5 dt 2 title "skam-gass-nonu-centring uniform", \
     'results/sknuc-nonuniform.dat' using 1:3 lc 5 dt 4 title "skam-gass-nonu-centring nonuniform", \
     x * 5e-1 lc rgbcolor "black" dt 1 lw 2 title '1st order', \
     x**2 * 1e1 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**3 * 1e2 lc rgbcolor "black" dt 3 lw 2 title '3rd order'
