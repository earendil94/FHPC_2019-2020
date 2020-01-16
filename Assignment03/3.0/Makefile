OB: openBody.c
	gcc -fopenmp -o openBody openBody.c -std=c99 -lrt -lm

IC: initialConditions.c
	mpicc -o initialConditions initialConditions.c -std=gnu99

MB: MPIBody.c
	mpicc -o MPIBody MPIBody.c -std=gnu99 -lm

DebugOB: openBody.c
	gcc -DDEBUG -fopenmp -o openBody openBody.c -std=c99 -lrt -lm

DebugNonOpen: openBody.c
	gcc -DDEBUG -o openBody openBody.c -std=c99 -lrt

DebugIC: initialConditions.c
	mpicc -DIC_DEBUG -o initialConditions initialConditions.c -std=gnu99

DebugMB: MPIBody.c
	mpicc -DDEBUG -o MPIBody MPIBody.c -std=gnu99 -lm

clean:
	rm *.ic 