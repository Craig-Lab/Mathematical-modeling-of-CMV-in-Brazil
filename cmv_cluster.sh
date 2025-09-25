#!/bin/bash
#SBATCH --account=def-craigm  # replace this with your own account
#SBATCH --ntasks=100             # number of processes
#SBATCH --mem-per-cpu=7G      # memory; default unit is megabytes
#SBATCH --time=5-00:00:00         # time (HH:MM:SS)

module load StdEnv/2020
module load r/4.2.2

# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
Rscript cedar_lenormand.r