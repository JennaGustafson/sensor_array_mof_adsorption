#!/bin/bash

output_dir=$1
regexp='Average loading absolute \[milligram\/gram framework\]'

grep $regexp $output_dir/*00.data | split_row.py
