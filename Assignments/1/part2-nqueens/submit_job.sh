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
printf "Setting stack size to unlimited..."
ulimit -s unlimited
ulimit -l unlimited
ulimit -a
printf

printf "Running Job...!"
printf "==============================================================================="
printf "Running compiled binary..."

# Module commands for compilation:
# module avail # shows all available modules in KUACC
# module list #list currently loaded modules.squeue
module load intel/ipsxe2019-u1ce # loads Intel compiler
# module load gcc/7.2.1/gcc # loads GNU compiler

#Enabling OMP settings
export OMP_NESTED=true
export OMP_MAX_ACTIVE_LEVELS=3
export OMP_NUM_THREADS=32

printf "\nExperiment 1 - Compact Affinity\n"
export KMP_AFFINITY=granularity=fine,compact # i removed "verbose"

printf "\n1.1 - Vanilla Parallelism"
printf "Serial, size 14"
./nqn_ser
printf "Parallel, size 14, 2 threads"
./nqn_par -t 2
printf "Parallel, size 14, 4 threads"
./nqn_par -t 4
printf "Parallel, size 14, 8 threads"
./nqn_par -t 8
printf "Parallel, size 14, 16 threads"
./nqn_par -t 16
printf "Parallel, size 14, 32 threads"
./nqn_par -t 32

printf "\n1.2 - Parallelism with Cutoff"
printf "Serial, size 14"
./nqn_ser
printf "Parallel, size 14, 2 threads, cutoff 5"
./nqn_par_c -t 2 -c 5
printf "Parallel, size 14, 4 threads, cutoff 5"
./nqn_par_c -t 4 -c 5
printf "Parallel, size 14, 8 threads, cutoff 5"
./nqn_par_c -t 8 -c 5
printf "Parallel, size 14, 16 threads, cutoff 5"
./nqn_par_c -t 16 -c 5
printf "Parallel, size 14, 32 threads, cutoff 5"
./nqn_par_c -t 32 -c 5

printf "\n1.3 - Early Stopping"
printf "Serial, size 14, early stop"
./nqn_ser_e
printf "Parallel, size 14, 2 threads, early stop"
./nqn_par_e -t 2
printf "Parallel, size 14, 4 threads, early stop"
./nqn_par_e -t 4
printf "Parallel, size 14, 8 threads, early stop"
./nqn_par_e -t 8
printf "Parallel, size 14, 16 threads, early stop"
./nqn_par_e -t 16
printf "Parallel, size 14, 32 threads, early stop"
./nqn_par_e -t 32

printf "\n1.4 - Problem Sizes"
printf "Parallel, size 13, 32 threads"
./nqn_par_13 -t 32
printf "Parallel, size 14, 32 threads"
./nqn_par -t 32
printf "Parallel, size 15, 32 threads"
./nqn_par_15 -t 32

printf "\nExperiment 2 - Scatter Affinity\n"
export KMP_AFFINITY=granularity=fine,scatter # i removed "verbose"

printf "\n2.1 - Vanilla Parallelism"
printf "Serial, size 14"
./nqn_ser
printf "Parallel, size 14, 2 threads"
./nqn_par -t 2
printf "Parallel, size 14, 4 threads"
./nqn_par -t 4
printf "Parallel, size 14, 8 threads"
./nqn_par -t 8
printf "Parallel, size 14, 16 threads"
./nqn_par -t 16
printf "Parallel, size 14, 32 threads"
./nqn_par -t 32

printf "\n2.2 - Parallelism with Cutoff"
printf "Serial, size 14"
./nqn_ser
printf "Parallel, size 14, 2 threads, cutoff 5"
./nqn_par_c -t 2 -c 5
printf "Parallel, size 14, 4 threads, cutoff 5"
./nqn_par_c -t 4 -c 5
printf "Parallel, size 14, 8 threads, cutoff 5"
./nqn_par_c -t 8 -c 5
printf "Parallel, size 14, 16 threads, cutoff 5"
./nqn_par_c -t 16 -c 5
printf "Parallel, size 14, 32 threads, cutoff 5"
./nqn_par_c -t 32 -c 5

printf "\n2.3 - Early Stopping"
printf "Serial, size 14, early stop"
./nqn_ser_e
printf "Parallel, size 14, 2 threads, early stop"
./nqn_par_e -t 2
printf "Parallel, size 14, 4 threads, early stop"
./nqn_par_e -t 4
printf "Parallel, size 14, 8 threads, early stop"
./nqn_par_e -t 8
printf "Parallel, size 14, 16 threads, early stop"
./nqn_par_e -t 16
printf "Parallel, size 14, 32 threads, early stop"
./nqn_par_e -t 32

printf "\n2.4 - Problem Sizes"
printf "Parallel, size 13, 32 threads"
./nqn_par_13 -t 32
printf "Parallel, size 14, 32 threads"
./nqn_par -t 32
printf "Parallel, size 15, 32 threads"
./nqn_par_15 -t 32