#!/bin/bash

output_dir=$1
line_regexp="Average loading absolute \[milligram\/gram framework\].*+\/-"
value_regexp="[0-9\.]*"

total_mass=0.0

for gas_mass in $(
  grep --only-matching "$line_regexp" $output_dir |
    grep --only-matching "$value_regexp" ); do

  total_mass=$( echo "$total_mass + $gas_mass" | bc -l )
done

echo $total_mass
