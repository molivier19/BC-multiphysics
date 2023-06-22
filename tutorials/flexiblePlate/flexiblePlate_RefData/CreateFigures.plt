# Gnuplot script


set style line 1 dt 1 lc rgb "red" lw 3
set style line 2 dt 2 lc rgb "red" lw 3
set style line 3 dt 3 lc rgb "red" lw 3
set style line 4 dt 1 lc rgb "blue" lw 3
set style line 5 dt 2 lc rgb "blue" lw 3
set style line 6 dt 3 lc rgb "blue" lw 3
set style line 7 dt 1 lc rgb "black" lw 3

#set datafile separator ";"

set xlabel "S/h coordinate"

set term qt 1

set ylabel "q/Q_s"
set key bottom left

stats "../postProcessing/singleGraph/solid/0.001/line_D_grad(T).xy" u (q_s0001=$5) every ::0::0 nooutput
stats "../postProcessing/singleGraph/solid/0.001/line_D_grad(T).xy" u (q_s001=$5) every ::0::0 nooutput
stats "../postProcessing/singleGraph/solid/0.01/line_D_grad(T).xy" u (q_s01=$5) every ::0::0 nooutput

plot \
    "RigidData_0001s.dat"    u 1:2 w l ls 1 title "Rigid plate, t = 0.0001s", \
    "RigidData_001s.dat"     u 1:2 w l ls 2 title "Rigid plate, t = 0.001s", \
    "RigidData_01s.dat"      u 1:2 w l ls 3 title "Rigid plate, t = 0.01s", \
    "FlexibleData_0001s.dat" u 1:2 w l ls 4 title "Flexible plate, t = 0.0001s", \
    "FlexibleData_001s.dat"  u 1:2 w l ls 5 title "Flexible plate, t = 0.001s", \
    "FlexibleData_01s.dat"   u 1:2 w l ls 6 title "Flexible plate, t = 0.01s", \
    "../postProcessing/singleGraph/solid/0.0001/line_D_grad(T).xy" u 1:($5/q_s0001) w l ls 7 title "Current simulation", \
    "../postProcessing/singleGraph/solid/0.001/line_D_grad(T).xy" u 1:($5/q_s001) w l ls 7 notitle, \
    "../postProcessing/singleGraph/solid/0.01/line_D_grad(T).xy" u 1:($5/q_s01) w l ls 7 notitle


set term qt 2

set ylabel "Normalize displacement magnitude (D/h)"
set key top left

plot \
    "RigidData_0001s.dat"    u 1:3 w l ls 1 title "Rigid plate, t = 0.0001s", \
    "RigidData_001s.dat"     u 1:3 w l ls 2 title "Rigid plate, t = 0.001s", \
    "RigidData_01s.dat"      u 1:3 w l ls 3 title "Rigid plate, t = 0.01s", \
    "FlexibleData_0001s.dat" u 1:3 w l ls 4 title "Flexible plate, t = 0.0001s", \
    "FlexibleData_001s.dat"  u 1:3 w l ls 5 title "Flexible plate, t = 0.001s", \
    "FlexibleData_01s.dat"   u 1:3 w l ls 6 title "Flexible plate, t = 0.01s", \
    "../postProcessing/singleGraph/solid/0.0001/line_D_grad(T).xy" u 1:(sqrt($2*$2+$3*$3)) w l ls 7 title "Current simulation", \
    "../postProcessing/singleGraph/solid/0.001/line_D_grad(T).xy"  u 1:(sqrt($2*$2+$3*$3)) w l ls 7 notitle, \
    "../postProcessing/singleGraph/solid/0.01/line_D_grad(T).xy"   u 1:(sqrt($2*$2+$3*$3)) w l ls 7 notitle

