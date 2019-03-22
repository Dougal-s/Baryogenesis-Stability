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
set output "time_ld=".string_ld.".png"
set multiplot title "l_d=".string_ld."; b:phi, y:l, g:h_u, c:h_d, b:sqrt(0.5*l_d)*d, r:0.5*log_{10}(T^2+1)"
set xtics out (0,100,200,300,400,500)
set ytics out (0,1,2)
set xlabel "time"
set ylabel "field"
set xrange [0:500]
set yrange [0:2]
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:(sqrt($3**2+$4**2)) lc rgb "#000000" pt 0
unset tics
set xlabel
set ylabel
plot "Phi.bin" binary format="%float32%uint32%26float32" record=0 using 1:(0.5*log10($13**2+$14**2+$15**2+$16**2+$17**2+$18**2+$19**2+$20**2+$21**2+$22**2+$23**2+$24**2+$25**2+$26**2+$27**2+$28**2+1)) lc rgb "#FF0000" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:(sqrt($5**2+$6**2)) lc rgb "#00FF00" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:(sqrt($9**2+$10**2)) lc rgb "#FFFF00" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:(sqrt($7**2+$8**2)) lc rgb "#00FFFF" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using 1:(sqrt(0.5*ld)*(sqrt($11**2+$12**2))) lc rgb "#0000FF" pt 0

