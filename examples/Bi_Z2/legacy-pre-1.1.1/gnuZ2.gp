set term post portrait  enhanced color "Helvetica,15"
set output 'Z2.eps'

set yrange [-.5:.5]
set xrange [-.5:.5]
set xtics -.5,.08333333333333333333,.5
set ytics -.5,.08333333333333333333,.5

set format y "%.1f";set format x "%.1f";set format z "%.1f"

set palette defined (-.9 'blue',0 'white', .9 'red' )

set arrow from -0.5,0 to 0.5,0 nohead lt 1 lw 3 lc rgb 'black'
set arrow from 0,-0.5 to 0,0.5 nohead lt 1 lw 3 lc rgb 'black'
set arrow from -0.5,.25 to 0.5,.25 nohead lt 1 lw 5 lc rgb 'black'
set arrow from 0.25,-0.5 to 0.25,0.5 nohead lt 1 lw 5 lc rgb 'black'
set arrow from -0.5,-.25 to 0.5,-.25 nohead lt 1 lw 5 lc rgb 'black'
set arrow from -0.25,-0.5 to -0.25,0.5 nohead lt 1 lw 5 lc rgb 'black'

set size square
set grid
set style fill transparent solid 0.5

plot 'NFIELD.dat' u 5:6:4 w p ps 2 lc palette pt 7 noti
