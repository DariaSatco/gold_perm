set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_perm_re.eps'
set autoscale 
set xtic auto
set xrange [0.2:5]                      
set ytic auto
set logscale y
set logscale x
set xlabel "{/Symbol w}, eV" font "Times-Roman, 20"
set ylabel "Re {/Symbol e} ({/Symbol w})" font "Times-Roman, 20" 
set key right top	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2 pt 7 ps 1.3
set style line 2 lt 1 lc rgb "#008040" lw 2 pt 11 ps 1.3
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 13 ps 1.3
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 2 ps 1
set style line 5 lt 1 lc rgb "#191970" lw 2 pt 7 ps 1
set style line 6 lt 1 lc rgb "#006400" lw 2 pt 64 ps 1
set style line 7 lt 0 lc rgb "#f055f0" lw 2 pt 65 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 2 pt 73 ps 1

plot "eps(w)_re.txt" u 1:(-($2)) smooth csplines w l ls 1 t "Drude-Lorentz",\
"eps(w)_re.txt" u 1:(-($3)) smooth csplines \
w l ls 3 t "Generalized Plasma",\
"eps(w)_re.txt" u 1:(-($4)) smooth csplines \
w l ls 4 t "Brendel-Bormann",\
"resultAu_eV.txt" u 1:(-($2)) w p ls 5 t "Palik data"
