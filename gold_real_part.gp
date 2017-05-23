set terminal postscript landscape color enhanced "Times-Roman" 16
set output 'gold_real_part_m.eps'
set autoscale 
set xtic auto
set xrange [1.0:35]    
set yrange [-1:10]                  
set ytic auto
set xtic (1,5,10,15,20,25,20,25,30,35)
set xlabel "{/Symbol w}, eV" font "Times-Roman, 20"
set ylabel "Re {/Symbol e} ({/Symbol w})" font "Times-Roman, 20" 
set key right top	#set the legend
# unset key

set style line 1 lt 1 lc rgb "#00eeee" lw 2 pt 7 ps 1
set style line 2 lt 1 lc rgb "#008040" lw 2 pt 11 ps 1
set style line 3 lt 1 lc rgb "#ffa040" lw 2 pt 13 ps 1
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 7 ps 0.8
set style line 5 lt 1 lc rgb "#191970" lw 2.5 pt 7 ps 1
set style line 6 lt 1 lc rgb "#006400" lw 2 pt 64 ps 1
set style line 7 lt 1 lc rgb "#ff4500" lw 2 pt 7 ps 0.8
set style line 8 lt 1 lc rgb "#ff1493" lw 2 pt 73 ps 1

plot "eps_minus_drude_re_im.txt" u 1:2 w p ls 4 t "Experimental data",\
"eps_minus_drude_re_im.txt" u 1:6 w l ls 8 t "Kramers-Kronig transform of Im {/Symbol e} ({/Symbol w})"

#"eps_minus_drude_re_im.txt" u 1:4 w l ls 5 t "Gauss approximation",\