# Gnuplot script

set style line 1 lt rgb "red" lw 4
set style line 2 lt rgb "blue" lw 2
set style line 3 lt rgb "black" lw 1
set style line 4 lt rgb "cyan" lw 1

# set datafile separator ";"

set xlabel "Theta (degrees)"
set ylabel "q/q_0"

time = "2"

heatMixed = sprintf("../mixed/postProcessing/sample_fluid/%s/wallHeatFluxKappa_InterfaceFluid.raw", time)
heatDN = sprintf("../DirichletNeumann/postProcessing/sample_fluid/%s/wallHeatFluxKappa_InterfaceFluid.raw", time)
TMixed = sprintf("../mixed/postProcessing/sample_solid/%s/T_InterfaceSolid.raw", time)
TDN = sprintf("../DirichletNeumann/postProcessing/sample_solid/%s/T_InterfaceSolid.raw", time)

stats heatMixed u 4
minMixed = STATS_min

stats heatDN u 4
minDN = STATS_min

set term qt 1

plot \
    "HeatTransferExperimental.csv"       u 1:4 title "Exp. Wieting", \
    "HeatTransferNumerical-Hongpeng.csv" u 1:2 w l title "Num. Hongpeng et al.", \
    "HeatTransfer-DN-cold.csv"           u 3:2 w l title "Ref. data (Dirichlet-Neumann)", \
    heatMixed u ((180/pi)*atan2($2,-$1)):($4/minMixed) title "Current, mixed", \
    heatDN u ((180/pi)*atan2($2,-$1)):($4/minDN) title "Current, D.-N."

set term qt 2

set ylabel "T (K)"

plot \
    "HeatTransferExperimental.csv" u 1:($2*5/9) title "Exp. Wieting", \
    TMixed u ((180/pi)*atan2($2,-$1)):4 title "Current, mixed", \
    TDN    u ((180/pi)*atan2($2,-$1)):4 title "Current, D.-N.", \

