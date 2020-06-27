CC = g++
CFLAGS  = -std=c++0x -g -Wall -O3

default: noise_remover

noise_remover: noise_remover.o
	$(CC) $(CFLAGS) noise_remover.o -lm -o noise_remover
	
noise_remover.o: noise_remover.cpp
	$(CC) $(CFLAGS) -c  noise_remover.cpp

clean:
	rm -rf *.o noise_remover
