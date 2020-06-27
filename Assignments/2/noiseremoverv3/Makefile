CC = nvcc
CFLAGS  = -g -O3

default: noise_remover_v3

noise_remover_v3: noise_remover_v3.o
	$(CC) $(CFLAGS) noise_remover_v3.o -lm -o noise_remover_v3
	
noise_remover_v3.o: noise_remover_v3.cu
	$(CC) $(CFLAGS) -c  noise_remover_v3.cu

clean:
	rm -rf *.o noise_remover_v3
