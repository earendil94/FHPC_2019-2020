all: nBodyMPI.c
	mpicc -o nBodyMPI nBodyMPI.c -lm

debug: nBodyMPI.c
	mpicc -DDEBUG -o nBodyMPI nBodyMPI.c -lm

clean:
	rm Debug* ic*