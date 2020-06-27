#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=noise_removal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short
#SBATCH --gres=gpu:tesla_k20m:1
#SBATCH # --gres=gpu:tesla_k40m:1 # with gpu specified
#SBATCH --time=00:10:00
#SBATCH --output=noise_removal.out
#SBATCH --exclusive

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

echo " "
echo "All tests are on coffee.pgm with -iter 100"
echo " "
echo "##################################################################################################"
echo "Serial"
./noiseremoverv0/noise_remover -i ./noiseremoverv0/coffee.pgm -o denoised_coffee0.png -iter 100

echo " "
echo "##################################################################################################"
echo "Parallel V1 - All Global Access"
./noiseremoverv1/noise_remover_v1 -i ./noiseremoverv1/coffee.pgm -o denoised_coffee1.png -iter 100

echo " "
echo "##################################################################################################"
echo "Parallel V2 - Data reuse with in-thread registers"
./noiseremoverv2/noise_remover_v2 -i ./noiseremoverv2/coffee.pgm -o denoised_coffee2.png -iter 100

echo " "
echo "##################################################################################################"
echo "Parallel V3 - Shared memory on Kernel 2"
./noiseremoverv3/noise_remover_v3 -i ./noiseremoverv3/coffee.pgm -o denoised_coffee3.png -iter 100

echo " "
echo "##################################################################################################"
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 1, CK2: 1, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 1 -ctwok 1 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 1, CK2: 2, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 1 -ctwok 2 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 1, CK2: 1, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 1 -ctwok 1 -usepinned 1
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 1, CK2: 2, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 1 -ctwok 2 -usepinned 1

echo "##################################################################################################"
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 2, CK2: 1, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 2 -ctwok 1 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 2, CK2: 2, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 2 -ctwok 2 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 2, CK2: 1, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 2 -ctwok 1 -usepinned 1
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 2, CK2: 2, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 2 -ctwok 2 -usepinned 1

echo "##################################################################################################"
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 3, CK2: 1, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 3 -ctwok 1 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 3, CK2: 2, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 3 -ctwok 2 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 3, CK2: 1, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 3 -ctwok 1 -usepinned 1
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 3, CK2: 2, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 3 -ctwok 2 -usepinned 1

echo "##################################################################################################"
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 4, CK2: 1, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 4 -ctwok 1 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 4, CK2: 2, No Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 4 -ctwok 2 -usepinned 0
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 4, CK2: 1, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 4 -ctwok 1 -usepinned 1
echo "Parallel V4 - Shared Memory and extra optimizations, CK1: 4, CK2: 2, Yes Pinned"
./noiseremoverv4/noise_remover_v4 -i ./noiseremoverv4/coffee.pgm -o denoised_coffee4.png -iter 100 -conek 4 -ctwok 2 -usepinned 1