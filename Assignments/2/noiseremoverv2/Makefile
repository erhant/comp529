CC = nvcc
CFLAGS  = -g -O3

default: noise_remover_v2

noise_remover_v2: noise_remover_v2.o
	$(CC) $(CFLAGS) noise_remover_v2.o -lm -o noise_remover_v2
	
noise_remover_v2.o: noise_remover_v2.cu
	$(CC) $(CFLAGS) -c  noise_remover_v2.cu

clean:
	rm -rf *.o noise_remover_v2
