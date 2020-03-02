#!bin/bash

# Create file to store Timing data in and save start time 
sudo touch Timing.txt
START=$(date)

# Run R script 
sudo Rscript 1_Simulation.R
END=$(date)

# Save start and end time in file
sudo echo "The start time was: ${START}" >> Timing.txt
sudo echo "The end time was: ${END}" >> Timing.txt

# If successful do git  pull and push
sudo git add .
sudo git commit -m "Simulation results"
sudo git pull
sudo git add .
sudo git commit -m "Manual merge"
sudo git push 
