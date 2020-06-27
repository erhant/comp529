#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=game-of-life-jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --time=01:50:00
#SBATCH --output=game-of-life-1.out

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
module load intel/ipsxe2019-u1ce # loads Intel compiler
module load gcc/7.2.1/gcc # loads GNU compiler

#Enabling OMP settings

#export KMP_AFFINITY=granularity=fine,compact # i removed "verbose"

#Serial version
echo "Running Serial"
## Remove log file
#rm log_s_1.txt
export OMP_NESTED=true
export OMP_MAX_ACTIVE_LEVELS=3
export OMP_NUM_THREADS=32
export KMP_AFFINITY=granularity=fine,compact # i removed "verbose"

./life -s 1 -n 2000 -i 500 -p 0.2 -d -l -t 1
# without the -l option, this will output a log_s_1234.txt file for the output of each iteration.
# s stands for "serial" and 1234 is the seed number.

#Parallel version
echo "Running first experiment: Strong Scaling"
echo "Fixed problem size (2000), threads 2,4,8,16,32"

./life -s 1 -n 2000 -i 500 -l -p 0.2 -d -t 2

./life -s 1 -n 2000 -i 500 -l -p 0.2 -d -t 4

./life -s 1 -n 2000 -i 500 -l -p 0.2 -d -t 8

./life -s 1 -n 2000 -i 500 -l -p 0.2 -d -t 16

# Remove log file
#rm log_p_1.txt
./life -s 1 -n 2000 -i 500 -l -p 0.2 -d -t 32

# Compare results
#echo "Running diff once"
#diff log_s_1.txt log_p_1.txt

#Parallel version
echo "Running second experiment: Problem Size"
echo "Fixed thread count (16), problem sizes 2000,4000,6000,8000,10000"
./life -s 2 -n 2000 -i 500 -p 0.2 -d -l -t 16
./life -s 2 -n 4000 -i 500 -p 0.2 -d -l -t 16
./life -s 2 -n 6000 -i 500 -p 0.2 -d -l -t 16
./life -s 2 -n 8000 -i 500 -p 0.2 -d -l -t 16
./life -s 2 -n 10000 -i 500 -p 0.2 -d -l -t 16