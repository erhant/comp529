#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=noise_coop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short
#SBATCH --gres=gpu:tesla_v100:1 # with gpu specified
#SBATCH --time=00:30:00
#SBATCH --output=noise_coop.out

### Some gpus and their compute capabilities:
### tesla_k20m -> 3.5 (?)
### tesla_k80  -> 3.7
### tesla_v100 -> 7.0
### tesla_k40m -> 3.5
### gtx_1080ti -> 6.1

################################################################################
##################### !!! DO NOT EDIT ABOVE THIS LINE !!! ######################
################################################################################
# Set stack size to unlimited
echo "Setting stack size to unlimited..."
ulimit -s unlimited
ulimit -l unlimited
ulimit -a
echo

echo "Running Job...!"
echo "==============================================================================="
echo "Running compiled binary..."

# Module commands for compilation:
# module avail # shows all available modules in KUACC
# module list #list currently loaded modules.squeue
module load cuda/10.1 # loads Intel compiler
module load gcc/7.2.1/gcc # loads GNU compiler

echo "Cooperative Groups Test"
./noise_remover_v3 -i coffee.pgm -o denoised_coffee.png