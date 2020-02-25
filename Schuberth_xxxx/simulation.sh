#!bin/bash

sudo Rscript 1_Simulation.R

git pull
git add .
git commit -m "Test simulation"
git push 
