/* Created by S. Cozzini and G.P. Brandino for the Course of Data Management @ MHPC
 * Last Revision by R. Aversa: October 2018
 */
#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#define FILESIZE 256 
int main(int argc, char **argv){

	int i,rank, size, bufsize, nints,offset;
	int *buf;

	MPI_File fh;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	bufsize = FILESIZE/size;
	nints = bufsize/sizeof(int);
	buf=(int*) malloc(nints*sizeof(int));

	offset=rank*bufsize;

	MPI_File_open(MPI_COMM_WORLD,"datafile",MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
	MPI_File_seek(fh, offset, MPI_SEEK_SET);
	MPI_File_read(fh, buf, nints, MPI_INT, &status);

	for (i=0;i<nints;i++){
		printf("rank: %d, buf[%d]: %d \n", rank, i, buf[i]);
	}
	free(buf);
	MPI_File_close(&fh);
	MPI_Finalize();
	
	return 0;
}
