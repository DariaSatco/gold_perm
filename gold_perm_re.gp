set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_perm_re1.eps'
set autoscale                 
set ytic auto
set xtic auto
set logscale y
set logscale x
set xlabel "{/Symbol w}, eV" font "Times-Roman, 20"
set ylabel "Re {/Symbol e} ({/Symbol w})" font "Times-Roman, 20" 
set key right top	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2.5 pt 7 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 1.5 pt 2 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 5 ps 1
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 5 ps 1
set style line 5 lt 1 lc rgb "#191970" lw 2 pt 7 ps 1
set style line 6 lt 1 lc rgb "#ff0000" lw 2 pt 2 ps 1.5
set style line 7 lt 0 lc rgb "#0080ff" lw 7 pt 7 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 7 pt 5 ps 1

plot "eps(w)_re.txt" u 1:(abs($5)) smooth csplines \
w l ls 3 t "Marachevsky",\
"eps(w)_re.txt" u 1:(abs($3)) smooth csplines \
w l ls 4 t "Generalized Plasma",\
"eps(w)_re.txt" u 1:(abs($2)) smooth csplines \
w l ls 5 t "Drude-Lorentz",\
"eps(w)_re.txt" u 1:(abs($4)) smooth csplines \
w l ls 6 t "Brendel-Bormann",\
"eps(w)_re.txt" u 1:(abs($6)) smooth csplines \
w l ls 8 t "Drude",\
"gold_eps_im_re_ev+Olmon.txt" u 1:(abs($2)) w p ls 7 t "Palik data"
