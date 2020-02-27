#!bin/bash

# Create file and save start time 
sudo touch Timing.txt
START=$(date)

# Run R script 
sudo Rscript 1_Simulation.R
END=$(date)

# Save start and end time in file

sudo echo "The start time was: ${START}" >> Timing.txt
sudo echo "The end time was: ${END}" >> Timing.txt

# If succesfull do git pull and push
git pull
git add .
git commit -m "Simulation results QWF3-1"
git push 
