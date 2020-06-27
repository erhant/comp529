#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=nqueens-job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=short
#SBATCH --time=01:50:00
#SBATCH --output=nqueens.out

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
# module load gcc/7.2.1/gcc # loads GNU compiler

#Enabling OMP settings
export OMP_NESTED=true
export OMP_MAX_ACTIVE_LEVELS=3
export OMP_NUM_THREADS=32

echo "Compact Affinity"
export KMP_AFFINITY=granularity=fine,compact # i removed "verbose"

echo "Serial, size 14"
./nqn_ser
echo "Parallel, size 14, 2 threads, cutoff 2"
./nqn_par_c -t 2 -c 2
echo "Parallel, size 14, 4 threads, cutoff 2"
./nqn_par_c -t 4 -c 2
echo "Parallel, size 14, 8 threads, cutoff 2"
./nqn_par_c -t 8 -c 2
echo "Parallel, size 14, 16 threads, cutoff 2"
./nqn_par_c -t 16 -c 2
echo "Parallel, size 14, 32 threads, cutoff 2"
./nqn_par_c -t 32 -c 2

echo "Parallel, size 14, 2 threads, cutoff 4"
./nqn_par_c -t 2 -c 4
echo "Parallel, size 14, 4 threads, cutoff 4"
./nqn_par_c -t 4 -c 4
echo "Parallel, size 14, 8 threads, cutoff 4"
./nqn_par_c -t 8 -c 4
echo "Parallel, size 14, 16 threads, cutoff 4"
./nqn_par_c -t 16 -c 4
echo "Parallel, size 14, 32 threads, cutoff 4"
./nqn_par_c -t 32 -c 4

echo "Parallel, size 14, 2 threads, cutoff 5"
./nqn_par_c -t 2 -c 5
echo "Parallel, size 14, 4 threads, cutoff 5"
./nqn_par_c -t 4 -c 5
echo "Parallel, size 14, 8 threads, cutoff 5"
./nqn_par_c -t 8 -c 5
echo "Parallel, size 14, 16 threads, cutoff 5"
./nqn_par_c -t 16 -c 5
echo "Parallel, size 14, 32 threads, cutoff 5"
./nqn_par_c -t 32 -c 5

echo "Parallel, size 14, 2 threads, cutoff 6"
./nqn_par_c -t 2 -c 6
echo "Parallel, size 14, 4 threads, cutoff 6"
./nqn_par_c -t 4 -c 6
echo "Parallel, size 14, 8 threads, cutoff 6"
./nqn_par_c -t 8 -c 6
echo "Parallel, size 14, 16 threads, cutoff 6"
./nqn_par_c -t 16 -c 6
echo "Parallel, size 14, 32 threads, cutoff 6"
./nqn_par_c -t 32 -c 6

echo "Parallel, size 14, 2 threads, cutoff 8"
./nqn_par_c -t 2 -c 8
echo "Parallel, size 14, 4 threads, cutoff 8"
./nqn_par_c -t 4 -c 8
echo "Parallel, size 14, 8 threads, cutoff 8"
./nqn_par_c -t 8 -c 8
echo "Parallel, size 14, 16 threads, cutoff 8"
./nqn_par_c -t 16 -c 8
echo "Parallel, size 14, 32 threads, cutoff 8"
./nqn_par_c -t 32 -c 8

echo "Parallel, size 14, 2 threads, cutoff 10"
./nqn_par_c -t 2 -c 10
echo "Parallel, size 14, 4 threads, cutoff 10"
./nqn_par_c -t 4 -c 10
echo "Parallel, size 14, 8 threads, cutoff 10"
./nqn_par_c -t 8 -c 10
echo "Parallel, size 14, 16 threads, cutoff 10"
./nqn_par_c -t 16 -c 10
echo "Parallel, size 14, 32 threads, cutoff 10"
./nqn_par_c -t 32 -c 10

echo "Scatter Affinity"
export KMP_AFFINITY=granularity=fine,scatter # i removed "verbose"

echo "Serial, size 14"
./nqn_ser
echo "Parallel, size 14, 2 threads, cutoff 2"
./nqn_par_c -t 2 -c 2
echo "Parallel, size 14, 4 threads, cutoff 2"
./nqn_par_c -t 4 -c 2
echo "Parallel, size 14, 8 threads, cutoff 2"
./nqn_par_c -t 8 -c 2
echo "Parallel, size 14, 16 threads, cutoff 2"
./nqn_par_c -t 16 -c 2
echo "Parallel, size 14, 32 threads, cutoff 2"
./nqn_par_c -t 32 -c 2

echo "Parallel, size 14, 2 threads, cutoff 4"
./nqn_par_c -t 2 -c 4
echo "Parallel, size 14, 4 threads, cutoff 4"
./nqn_par_c -t 4 -c 4
echo "Parallel, size 14, 8 threads, cutoff 4"
./nqn_par_c -t 8 -c 4
echo "Parallel, size 14, 16 threads, cutoff 4"
./nqn_par_c -t 16 -c 4
echo "Parallel, size 14, 32 threads, cutoff 4"
./nqn_par_c -t 32 -c 4

echo "Parallel, size 14, 2 threads, cutoff 5"
./nqn_par_c -t 2 -c 5
echo "Parallel, size 14, 4 threads, cutoff 5"
./nqn_par_c -t 4 -c 5
echo "Parallel, size 14, 8 threads, cutoff 5"
./nqn_par_c -t 8 -c 5
echo "Parallel, size 14, 16 threads, cutoff 5"
./nqn_par_c -t 16 -c 5
echo "Parallel, size 14, 32 threads, cutoff 5"
./nqn_par_c -t 32 -c 5

echo "Parallel, size 14, 2 threads, cutoff 6"
./nqn_par_c -t 2 -c 6
echo "Parallel, size 14, 4 threads, cutoff 6"
./nqn_par_c -t 4 -c 6
echo "Parallel, size 14, 8 threads, cutoff 6"
./nqn_par_c -t 8 -c 6
echo "Parallel, size 14, 16 threads, cutoff 6"
./nqn_par_c -t 16 -c 6
echo "Parallel, size 14, 32 threads, cutoff 6"
./nqn_par_c -t 32 -c 6

echo "Parallel, size 14, 2 threads, cutoff 8"
./nqn_par_c -t 2 -c 8
echo "Parallel, size 14, 4 threads, cutoff 8"
./nqn_par_c -t 4 -c 8
echo "Parallel, size 14, 8 threads, cutoff 8"
./nqn_par_c -t 8 -c 8
echo "Parallel, size 14, 16 threads, cutoff 8"
./nqn_par_c -t 16 -c 8
echo "Parallel, size 14, 32 threads, cutoff 8"
./nqn_par_c -t 32 -c 8

echo "Parallel, size 14, 2 threads, cutoff 10"
./nqn_par_c -t 2 -c 10
echo "Parallel, size 14, 4 threads, cutoff 10"
./nqn_par_c -t 4 -c 10
echo "Parallel, size 14, 8 threads, cutoff 10"
./nqn_par_c -t 8 -c 10
echo "Parallel, size 14, 16 threads, cutoff 10"
./nqn_par_c -t 16 -c 10
echo "Parallel, size 14, 32 threads, cutoff 10"
./nqn_par_c -t 32 -c 10