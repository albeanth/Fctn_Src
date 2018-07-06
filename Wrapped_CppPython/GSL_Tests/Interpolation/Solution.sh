reset
set output "GSL_Interpolate_Test.eps"
set term postscript eps enhanced color blacktext "Helvetica" 18
set grid ytics mytics
set key
f(x) = sin(x)
plot  [0:pi] f(x),\
"data.out" index 0 using 1:2 with points pt 1 title 'Spline'

reset
set output "GSL_Interpolate_Deriv_Test.eps"
set term postscript eps enhanced color blacktext "Helvetica" 18
set grid ytics mytics
set key
f(x) = cos(x)
plot  [0:pi] f(x),\
"data.out" index 0 using 1:3 with points pt 1 title 'Deriv Spline'
