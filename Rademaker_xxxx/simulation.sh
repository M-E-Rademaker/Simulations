#!bin/bash

# Create file and directory and save start time
sudo mkdir -p Data_simulation
sudo touch Timing.txt
START=$(date)

# Run R script
sudo Rscript 1_Simulation.R
END=$(date)

# Save start and end time in file

sudo echo "The start time was: ${START}" >> Timing.txt
sudo echo "The end time was: ${END}" >> Timing.txt

# If succesful, add, commit and  do git pull and push
sudo git add .
sudo git commit -m "Simulation_results_hpc"
sudo git pull
sudo git push
