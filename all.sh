#!/bin/bash

file=$1
t=$2

cat ${file} | timeout ${t} build/main methods/dopri5.txt 1 explicit
mv file.log explicit.log

cat ${file} | timeout ${t} build/main methods/dopri5.txt 1 implicit
mv file.log implicit.log

cat ${file} | timeout ${t} build/main methods/dopri5.txt 1 adams
mv file.log adams.log

cat ${file} | timeout ${t} build/main methods/rk4.txt 1 adaptive
mv file.log adaptive.log

python3 plot.py
mv trajectory_comparison.png "plots/$(basename ${file}).png"
