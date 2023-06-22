# Gnuplot script

set style line 1 lw 4 dt 4 lc rgb "#00FF0000"   
set style line 2 lw 2 dt 3 lc rgb "#00CC0000"   
set style line 3 lw 1 dt 2 lc rgb "#00990000"   
set style line 4 lw 1 dt 1 lc rgb "#00660000"   
set style line 5 lw 2 dt 1 lc rgb "blue"  
set style line 6 lw 1 dt 1 lc rgb "black" 

set xlabel "Time (s)"
set ylabel "Net volumetric influx (m^3/s)"
set xrange [0:0.01]

plot \
    "Degroote_Dimentional.csv"           u 1:2 w l ls 6 title "Ref. Degroote et al.", \
    "InletOutletFlowRate-veryCoarse.dat" u 1:(-$2-$3) w l ls 1 title "Ref. very coarse mesh", \
    "InletOutletFlowRate-coarse.dat"     u 1:(-$2-$3) w l ls 2 title "Ref. coarse mesh", \
    "InletOutletFlowRate-medium.dat"     u 1:(-$2-$3) w l ls 3 title "Ref. medium mesh", \
    "InletOutletFlowRate-fine.dat"       u 1:(-$2-$3) w l ls 4 title "Ref. fine mesh", \
    "../../Step/InletOutletFlowRate.dat" u 1:(-$2-$3) w l ls 5 title "Current 'Step' case"

