set style data linespoints

set logscale
set format x "$10^{%L}$"
set format y "$10^{%L}$"

set xlabel "mean dx"
set ylabel "error"

set xrange [1e-1:1e-3]

unset key

plot 'results/uniform.dat' using 1:2 lc 1 dt 1 title "uniform l2" at end, \
     'results/uniform.dat' using 1:3 lc 1 dt 2 title "uniform linf" at end, \
     'results/nonuniform.dat' using 1:2 lc 2 dt 1 title "nonuniform l2" at end, \
     'results/nonuniform.dat' using 1:3 lc 2 dt 2 title "nonuniform linf" at end, \
     x * 5e-1 lc rgbcolor "black" dt 1 lw 2 title '1st order', \
     x**2 * 1e1 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**3 * 1e2 lc rgbcolor "black" dt 3 lw 2 title '3rd order'
