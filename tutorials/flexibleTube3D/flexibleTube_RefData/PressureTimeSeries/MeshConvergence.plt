# Gnuplot script

set style line 1 lw 4 dt 3 lc rgb "#00FF0000"   
set style line 2 lw 2 dt 2 lc rgb "#00CC0000"   
set style line 3 lw 1 dt 1 lc rgb "#00990000"   
set style line 4 lw 1 dt 1 lc rgb "#00660000"   
set style line 5 lw 2 dt 1 lc rgb "blue"  
set style line 6 lw 1 dt 1 lc rgb "black" 

set xlabel "Time (s)"
set ylabel "p/p_{pulse}"
set xrange [0:0.01]
# set datafile separator ";"

plot \
    "Pressure_VeryCoarse.txt" u 1:2 w l ls 1 title "Ref. very coarse mesh", \
    "Pressure_Coarse.txt"     u 1:2 w l ls 2 title "Ref. coarse mesh", \
    "Pressure_Medium.txt"     u 1:2 w l ls 3 title "Ref. medium mesh", \
    "Pressure_Fine.txt"       u 1:2 w l ls 4 title "Ref. fine mesh", \
    "Pressure_VeryCoarse.txt" u 1:3 w l ls 1 notitle, \
    "Pressure_Coarse.txt"     u 1:3 w l ls 2 notitle, \
    "Pressure_Medium.txt"     u 1:3 w l ls 3 notitle, \
    "Pressure_Fine.txt"       u 1:3 w l ls 4 notitle, \
    "Pressure_VeryCoarse.txt" u 1:4 w l ls 1 notitle, \
    "Pressure_Coarse.txt"     u 1:4 w l ls 2 notitle, \
    "Pressure_Medium.txt"     u 1:4 w l ls 3 notitle, \
    "Pressure_Fine.txt"       u 1:4 w l ls 4 notitle, \
    "../../Smooth/postProcessing/PressureProbes/fluid/0/p" u 1:2 ls 5 w l title "Current 'Smooth' case", \
    "../../Smooth/postProcessing/PressureProbes/fluid/0/p" u 1:3 ls 5 w l notitle, \
    "../../Smooth/postProcessing/PressureProbes/fluid/0/p" u 1:4 ls 5 w l notitle
