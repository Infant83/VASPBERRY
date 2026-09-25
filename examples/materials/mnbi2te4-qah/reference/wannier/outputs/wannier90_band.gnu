set style data dots
set nokey
set xrange [0: 2.28570]
set yrange [ -3.79049 :  8.47027]
set arrow from  0.83662,  -3.79049 to  0.83662,   8.47027 nohead
set arrow from  1.31965,  -3.79049 to  1.31965,   8.47027 nohead
set xtics ("G"  0.00000,"M"  0.83662,"K"  1.31965,"G"  2.28570)
 plot "wannier90_band.dat"
