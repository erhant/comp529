# Using intel compiler
icc -Wall -O3 -qopenmp -o nqn_ser nqueens_serial.c
icc -Wall -O3 -qopenmp -o nqn_ser_e nqueens_serial-early.c
icc -Wall -O3 -qopenmp -o nqn_par nqueens_parallel.c
icc -Wall -O3 -qopenmp -o nqn_par_13 nqueens_parallel_13.c
icc -Wall -O3 -qopenmp -o nqn_par_15 nqueens_parallel_15.c
icc -Wall -O3 -qopenmp -o nqn_par_e nqueens_parallel-early.c
icc -Wall -O3 -qopenmp -o nqn_par_c nqueens_parallel-cutoff.c