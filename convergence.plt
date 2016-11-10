set style data linespoints

set logscale
set format x "$10^{%L}$"
set format y "$10^{%L}$"

set xlabel "mean dx"
set ylabel "error"

set xrange [1e-1:1e-3]
set yrange [1e-7:1]

set key outside

plot 'results/sk-uniform.dat' using 1:2 lc 1 dt 1 title "sk uniform l2", \
     'results/sk-uniform.dat' using 1:3 lc 1 dt 2 title "sk uniform linf", \
     'results/sk-nonuniform.dat' using 1:2 lc 2 dt 1 title "sk nonuniform l2", \
     'results/sk-nonuniform.dat' using 1:3 lc 2 dt 2 title "sk nonuniform linf", \
     'results/ls-uniform.dat' using 1:2 lc 3 dt 1 title "ls uniform l2", \
     'results/ls-uniform.dat' using 1:3 lc 3 dt 2 title "ls uniform linf", \
     'results/ls-nonuniform.dat' using 1:2 lc 4 dt 1 title "ls nonuniform l2", \
     'results/ls-nonuniform.dat' using 1:3 lc 4 dt 2 title "ls nonuniform linf", \
     'results/cf-uniform.dat' using 1:2 lc 5 dt 1 title "cf uniform l2", \
     'results/cf-uniform.dat' using 1:3 lc 5 dt 2 title "cf uniform linf", \
     'results/cf-nonuniform.dat' using 1:2 lc 6 dt 1 title "cf nonuniform l2", \
     'results/cf-nonuniform.dat' using 1:3 lc 6 dt 2 title "cf nonuniform linf", \
     x * 5e-1 lc rgbcolor "black" dt 1 lw 2 title '1st order', \
     x**2 * 1e1 lc rgbcolor "black" dt 2 lw 2 title '2nd order', \
     x**3 * 1e2 lc rgbcolor "black" dt 3 lw 2 title '3rd order'
