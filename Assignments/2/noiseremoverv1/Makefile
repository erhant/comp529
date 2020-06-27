CC = nvcc
CFLAGS  = -g -O3

default: noise_remover_v1

noise_remover_v1: noise_remover_v1.o
	$(CC) $(CFLAGS) noise_remover_v1.o -lm -o noise_remover_v1
	
noise_remover_v1.o: noise_remover_v1.cu
	$(CC) $(CFLAGS) -c  noise_remover_v1.cu

clean:
	rm -rf *.o noise_remover_v1
