#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
# This is a performance study test for the 3rd assignment of Parallel Programming.
#
# -= Resources =-
#
#SBATCH --job-name=spmv-jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --output=spmv-job.out

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

module load gcc/9.1.0
module load mpich/3.2

# Print CPU Info
lscpu

echo " "
echo "##################################################################################################"

# Serial (version 0)
echo " "
echo "Serial version on Cube Coup dt16, 20 timesteps"
./V0/build/spmv Cube_Coup_dt6.mtx 20

# MPI (version 1)
echo " "
echo "V1 - MPI ONLY with 1 process on Cube Coup dt16, 20 timesteps"
mpiexec -n 1 ./V1/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V1 - MPI ONLY with 2 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 2 ./V1/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V1 - MPI ONLY with 4 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 4 ./V1/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V1 - MPI ONLY with 8 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 8 ./V1/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V1 - MPI ONLY with 16 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 16 ./V1/build/spmv Cube_Coup_dt6.mtx 20

# MPI + OpenMP (version 2)
export KMP_AFFINITY=verbose,granularity=fine,compact
export OMP_NUM_THREADS=16
echo " "
echo "V2 - MPI + OpenMP with 1 process on Cube Coup dt16, 20 timesteps"
mpiexec -n 1 ./V2/build/spmv Cube_Coup_dt6.mtx 20

export OMP_NUM_THREADS=8
echo " "
echo "V2 - MPI + OpenMP with 2 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 2 ./V2/build/spmv Cube_Coup_dt6.mtx 20

export OMP_NUM_THREADS=4
echo " "
echo "V2 - MPI + OpenMP with 4 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 4 ./V2/build/spmv Cube_Coup_dt6.mtx 20

export OMP_NUM_THREADS=2
echo " "
echo "V2 - MPI + OpenMP with 8 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 8 ./V2/build/spmv Cube_Coup_dt6.mtx 20

export OMP_NUM_THREADS=1
echo " "
echo "V2 - MPI + OpenMP with 16 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 16 ./V2/build/spmv Cube_Coup_dt6.mtx 20

# MPI with Load Balancing (version 3)
echo " "
echo "V3 - MPI with Load Balancing with 1 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 1 ./V3/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 2 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 2 ./V3/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 4 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 4 ./V3/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 8 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 8 ./V3/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 16 processes on Cube Coup dt16, 20 timesteps"
mpiexec -n 16 ./V3/build/spmv Cube_Coup_dt6.mtx 20

echo " "
echo "##################################################################################################"

# Serial (version 0)
echo " "
echo "Serial version on Flan 1565, 20 timesteps"
./V0/build/spmv Flan_1565.mtx 20

# MPI (version 1)
echo " "
echo "V1 - MPI ONLY with 1 process on Flan 1565, 20 timesteps"
mpiexec -n 1 ./V1/build/spmv Flan_1565.mtx 20

echo " "
echo "V1 - MPI ONLY with 2 processes on Flan 1565, 20 timesteps"
mpiexec -n 2 ./V1/build/spmv Flan_1565.mtx 20

echo " "
echo "V1 - MPI ONLY with 4 processes on Flan 1565, 20 timesteps"
mpiexec -n 4 ./V1/build/spmv Flan_1565.mtx 20

echo " "
echo "V1 - MPI ONLY with 8 processes on Flan 1565, 20 timesteps"
mpiexec -n 8 ./V1/build/spmv Flan_1565.mtx 20

echo " "
echo "V1 - MPI ONLY with 16 processes on Flan 1565, 20 timesteps"
mpiexec -n 16 ./V1/build/spmv Flan_1565.mtx 20

# MPI + OpenMP (version 2)
export KMP_AFFINITY=verbose,granularity=fine,compact
export OMP_NUM_THREADS=16
echo " "
echo "V2 - MPI + OpenMP with 1 process on Flan 1565, 20 timesteps"
mpiexec -n 1 ./V2/build/spmv Flan_1565.mtx 20

export OMP_NUM_THREADS=8
echo " "
echo "V2 - MPI + OpenMP with 2 processes on Flan 1565, 20 timesteps"
mpiexec -n 2 ./V2/build/spmv Flan_1565.mtx 20

export OMP_NUM_THREADS=4
echo " "
echo "V2 - MPI + OpenMP with 4 processes on Flan 1565, 20 timesteps"
mpiexec -n 4 ./V2/build/spmv Flan_1565.mtx 20

export OMP_NUM_THREADS=2
echo " "
echo "V2 - MPI + OpenMP with 8 processes on Flan 1565, 20 timesteps"
mpiexec -n 8 ./V2/build/spmv Flan_1565.mtx 20

export OMP_NUM_THREADS=1
echo " "
echo "V2 - MPI + OpenMP with 16 processes on Flan 1565, 20 timesteps"
mpiexec -n 16 ./V2/build/spmv Flan_1565.mtx 20

# MPI with Load Balancing (version 3)
echo " "
echo "V3 - MPI with Load Balancing with 1 process on Flan 1565, 20 timesteps"
mpiexec -n 1 ./V3/build/spmv Flan_1565.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 2 processes on Flan 1565, 20 timesteps"
mpiexec -n 2 ./V3/build/spmv Flan_1565.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 4 processes on Flan 1565, 20 timesteps"
mpiexec -n 4 ./V3/build/spmv Flan_1565.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 8 processes on Flan 1565, 20 timesteps"
mpiexec -n 8 ./V3/build/spmv Flan_1565.mtx 20

echo " "
echo "V3 - MPI with Load Balancing with 16 processes on Flan 1565, 20 timesteps"
mpiexec -n 16 ./V3/build/spmv Flan_1565.mtx 20