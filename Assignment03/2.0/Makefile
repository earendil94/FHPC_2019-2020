all: initialConditions.c
	mpicc -o initialConditions initialConditions.c

debug: initialConditions.c
	mpicc -DDEBUG -o initialConditions initialConditions.c

clean:
	rm *.x