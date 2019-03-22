set terminal png size 1600,1200
set key noautotitle
set title ''
set tmargin 2
set lmargin 8
set rmargin 3
set bmargin 4
unset multiplot

ld=1000 # put value of lambda_d here

set output "initial.png"
set multiplot title "t=0"
set xtics out (0,8192)
set ytics out (0,1,2)
set xlabel "lattice"
set ylabel "field"
set xrange [0:8192] # put lattice size here
set yrange [0:2]
plot "initialPhi.bin" binary format="%float32%uint32%26float32" record=0 using ($2):($13**2+$14**2+$15**2+$16**2+$17**2+$18**2+$19**2+$20**2+$21**2+$22**2+$23**2+$24**2+$25**2+$26**2+$27**2+$28**2) lc rgb "#FF0000" pt 0
unset tics
set xlabel
set ylabel
plot "initialPhi.bin" binary format="%float32%uint32%26float32" using ($2):(1E3*(sqrt($11**2+$12**2))) lc rgb "#0000FF" pt 0
plot "initialPhi.bin" binary format="%float32%uint32%26float32" using ($2):(sqrt($5**2+$6**2)) lc rgb "#00FF00" pt 0
plot "initialPhi.bin" binary format="%float32%uint32%26float32" using ($2):(sqrt($9**2+$10**2)) lc rgb "#FFFF00" pt 0
plot "initialPhi.bin" binary format="%float32%uint32%26float32" using ($2):(1E3*sqrt($7**2+$8**2)) lc rgb "#00FFFF" pt 0
plot "initialPhi.bin" binary format="%float32%uint32%26float32" using ($2):(1E5*sqrt($3**2+$4**2)) lc rgb "#000000" pt 0
unset multiplot
}
