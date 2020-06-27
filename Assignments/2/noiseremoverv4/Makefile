CC = nvcc
CFLAGS  = -g -O3

default: noise_remover_v4

noise_remover_v4: noise_remover_v4.o
	$(CC) $(CFLAGS) noise_remover_v4.o -lm -o noise_remover_v4
	
noise_remover_v4.o: noise_remover_v4.cu
	$(CC) $(CFLAGS) -c  noise_remover_v4.cu

clean:
	rm -rf *.o noise_remover_v4
