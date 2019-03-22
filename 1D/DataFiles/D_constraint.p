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

set output "D_constraint_ld=".string_ld.".png"
set multiplot title "l_d=".string_ld."; black:D, y:|l|^2, g:|h_u|^2, c:|h_d|^2, b:|d|^2"
set xtics out (0,100,200,300,400,500)
set ytics out (-1E-6,0,1E-6,2E-6)
set xlabel "time"
set ylabel "field"
set xrange [0:500]
set yrange [-1E-6:2E-6]
plot "Phi.bin" binary format="%float32%uint32%26float32" record = 0 using 1:($11**2+$12**2) lc rgb "#0000FF" pt 0
unset tics
set xlabel
set ylabel
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:($5**2+$6**2) lc rgb "#00FF00" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:($9**2+$10**2) lc rgb "#FFFF00" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:($7**2+$8**2) lc rgb "#00FFFF" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:($5**2+$6**2-$7**2-$8**2-$9**2-$10**2+0.5*$11**2+0.5*$12**2) lc rgb "#000000" pt 1
