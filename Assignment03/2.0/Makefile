all: initialConditions.c
	mpicc -o initialConditions initialConditions.c -std=gnu99

debug: initialConditions.c
	mpicc -DDEBUG -o initialConditions initialConditions.c -std=gnu99

clean:
	rm *.x *.kappa