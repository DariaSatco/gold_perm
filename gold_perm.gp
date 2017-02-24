set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_perm.eps'
set autoscale 
set xtic auto                      
set ytic auto
set logscale y
set xlabel "{/Symbol w}, eV" font "Times-Roman, 20"
set ylabel "{/Symbol e} (i{/Symbol w})" font "Times-Roman, 20" 
set key right top	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2 pt 7 ps 1.3
set style line 2 lt 1 lc rgb "#008040" lw 2 pt 11 ps 1.3
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 13 ps 1.3
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 2 ps 1
set style line 5 lt 1 lc rgb "#191970" lw 2 pt 1 ps 1
set style line 6 lt 1 lc rgb "#006400" lw 2 pt 64 ps 1
set style line 7 lt 0 lc rgb "#f055f0" lw 2 pt 65 ps 1
set style line 8 lt 0 lc rgb "#804014" lw 2 pt 73 ps 1

plot "golddata.txt" u 1:2 smooth csplines w l ls 1 t "Kramers-Kronig",\
"golddata.txt" u 1:4 smooth csplines \
w l ls 3 t "Marachevsky",\
"golddata.txt" u 1:5 smooth csplines \
w l ls 4 t "Generalized Plasma",\
"golddata.txt" u 1:6 smooth csplines \
w l ls 5 t "Brendel-Bormann",\
"golddata.txt" u 1:7 smooth csplines \
w l ls 6 t "Lorentz-Drude"
