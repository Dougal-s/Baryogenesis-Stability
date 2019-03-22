set terminal png size 1600,1200
set key noautotitle
set title ''
set tmargin 2
set lmargin 8
set rmargin 3
set bmargin 4
unset multiplot

ld=1000 # put value of lambda_d here

do for [time in "0 50 60 70 80 90 100 200 300 500"] { # choose time samples here
string_ld = gprintf("%g",ld)
string_time = gprintf("%g",time)
set output "space_ld=".string_ld."_t=".string_time.".png"
set multiplot title "l_d=".string_ld.", t=".string_time."
set xtics out (0,8192)
set ytics out (0,1,2)
set xlabel "lattice"
set ylabel "field"
set xrange [0:8192] # put lattice size here
set yrange [0:2]
plot "Phi.bin" binary format="%float32%uint32%26float32" record=0 using (($1==time)?($2):NaN):(0.5*log10($13**2+$14**2+$15**2+$16**2+$17**2+$18**2+$19**2+$20**2+$21**2+$22**2+$23**2+$24**2+$25**2+$26**2+$27**2+$28**2+1)) lc rgb "#FF0000" pt 0
unset tics
set xlabel
set ylabel
plot "Phi.bin" binary format="%float32%uint32%26float32" using (($1==time)?($2):NaN):(sqrt(0.5*ld)*(sqrt($11**2+$12**2))) lc rgb "#0000FF" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using (($1==time)?($2):NaN):(sqrt($5**2+$6**2)) lc rgb "#00FF00" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using (($1==time)?($2):NaN):(sqrt($9**2+$10**2)) lc rgb "#FFFF00" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using (($1==time)?($2):NaN):(sqrt($7**2+$8**2)) lc rgb "#00FFFF" pt 0
plot "Phi.bin" binary format="%float32%uint32%26float32" using (($1==time)?($2):NaN):(sqrt($3**2+$4**2)) lc rgb "#000000" pt 0
unset multiplot
}
