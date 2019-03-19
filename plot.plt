set term postscript eps enhanced color font "Arial, 24"


set palette defined (0 "white", 1 "red")
#set palette defined (0 "white", 0.5 "blue", 1 "red")
set cbrange [0:1]

set yrange [-1.3:1.3]
set xtics 0.2

set xlabel "Field strength (a.u.)"
set ylabel "Floquet quasienergy (a.u.)"
unset colorbox

unset key
set output "floquet_qeps.eps"
p for [i=1:2*(2*32+1)] "eps_Floquet.out" u 1:(column(i*2)):(column(i*2+1)) w l lt 1 lc palette lw 6




unset output