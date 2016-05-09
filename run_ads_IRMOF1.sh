#!/bin/bash

CO2_COMP=$1
CH4_COMP=$2
N2_COMP=$3
C2H6_COMP=$4
MOF=$5

mkdir -p outputs
mkdir -p inputs
cp simulation_template.file simulation.input
sed -i "s/^FrameworkName .*$/FrameworkName $MOF/" simulation.input
sed -i "20 s/[0-9]*\.[0-9]*/$CO2_COMP/" simulation.input
sed -i "31 s/[0-9]*\.[0-9]*/$CH4_COMP/" simulation.input
sed -i "42 s/[0-9]*\.[0-9]*/$N2_COMP/" simulation.input
sed -i "53 s/[0-9]*\.[0-9]*/$C2H6_COMP/" simulation.input

cp simulation.input sim_${CO2_COMP}_${CH4_COMP}_${N2_COMP}_${C2H6_COMP}.input
mv sim_${CO2_COMP}_${CH4_COMP}_${N2_COMP}_${C2H6_COMP}.input inputs
simulate simulation.input

mass_p1=$(./calculate_mass.sh ./Output/System_0/*00.data)
mass_p2=$(./calculate_mass.sh ./Output/System_0/*06.data)

mkdir Out_${CO2_COMP}_${CH4_COMP}_${N2_COMP}_${C2H6_COMP}
cp -r Output Out_${CO2_COMP}_${CH4_COMP}_${N2_COMP}_${C2H6_COMP}
mv Out_${CO2_COMP}_${CH4_COMP}_${N2_COMP}_${C2H6_COMP} outputs/

echo $mass_p1 $mass_p2
