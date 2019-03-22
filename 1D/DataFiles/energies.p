set terminal png size 800,600
set key noautotitle
set title ''
set tmargin 2
set lmargin 8
set rmargin 3
set bmargin 4
unset multiplot

ld=1000 # put value of lambda_d here

string_ld = gprintf("%g",ld)

set output "energies_ld=".string_ld.".png"
set multiplot title "l_d=".string_ld."; black:E, r:K, g:G, b:V"
set xtics out (0,100,200,300,400,500,600,700,800,900,1000)
set ytics out (-1E-2,0,1E-2)
set xlabel "time"
set ylabel "field"
set xrange [0:500]
set yrange [-5E-2:5E-2]
plot "energies.txt" using 1:2 lc rgb "#000000" pt 0
unset tics
set xlabel
set ylabel
plot "energies.txt" using 1:3 lc rgb "#FF0000" pt 0
plot "energies.txt" using 1:4 lc rgb "#00FF00" pt 0
plot "energies.txt" using 1:5 lc rgb "#0000FF" pt 0
unset multiplot

set output "energy_flow_ld=".string_ld.".png"
set multiplot title "l_d=".string_ld."; grey: E, black:E_{phi}, g:(K+G)_{AD}, b:(K+G)_d, r:(K+G)_{gluon}, y:V_{AD+d+gluon}"
set xtics out (0,100,200,300,400,500,600,700,800,900,1000)
set ytics out (-1E-3,-1E-4,0,1E-4,1E-3)
set xlabel "time"
set ylabel "field"
set xrange [0:500]
set yrange [-1E-3:1E-3]
plot "energies.txt" using 1:($6+$19+$32) lc rgb "#000000" pt 0
unset tics
set xlabel
set ylabel
plot "energies.txt" using 1:($7+$8+$9+$20+$21+$22) lc rgb "#00FF00" pt 0
plot "energies.txt" using 1:($10+$23) lc rgb "#0000FF" pt 0
plot "energies.txt" using 1:($11+$12+$13+$14+$15+$16+$17+$18+$24+$25+$26+$27+$28+$29+$30+$31) lc rgb "#FF0000" pt 0
plot "energies.txt" using 1:33 lc rgb "#FFFF00" pt 0
plot "energies.txt" using 1:2 lc rgb "#AAAAAA" pt 0
unset multiplot

