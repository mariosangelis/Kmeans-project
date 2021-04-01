/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_main.c   (an sequential version)                      */
/*   Description:  This program shows an example on how to call a subroutine */
/*                 that implements a simple k-means clustering algorithm     */
/*                 based on Euclid distance.                                 */
/*   Input file format:                                                      */
/*                 ascii  file: each line contains 1 data object             */
/*                 binary file: first 4-byte integer is the number of data   */
/*                 objects and 2nd integer is the no. of features (or        */
/*                 coordinates) of each object                               */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strtok() *
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     /* getopt() */

int      _debug;
#include "kmeans.h"

/*---< usage() >------------------------------------------------------------*/
static void usage(char *argv0, double threshold) {
char *help =
	"Usage: %s [switches] -i filename -n num_clusters\n"
	"       -i filename    : file containing data to be clustered\n"
	"       -c centers     : file containing initial centers. default: filename\n"
	"       -b             : input file is in binary format (default no)\n"
	"       -n num_clusters: number of clusters (K must > 1)\n"
	"       -t threshold   : threshold value (default %.4f)\n"
	"       -o             : output timing results (default no)\n"
	"       -q             : quiet mode\n"
	"       -d             : enable debug mode\n"
	"       -h             : print this help information\n"
	"       -f             : output file id(default 1)\n";
fprintf(stderr, help, argv0, threshold);
exit(-1);
}

/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
			int     opt;
			FILE *fd;
	extern char   *optarg;
	extern int     optind;
			int     i, j, isBinaryFile, is_output_timing, verbose;
			int     output_id;
			int     numClusters, numCoords, numObjs;
			int    *membership;    /* [numObjs] */
			char   *filename, *center_filename;
			double **objects;       /* [numObjs][numCoords] data objects */
			double **clusters;      /* [numClusters][numCoords] cluster center */
			float   threshold;
			double  timing, io_timing, clustering_timing;
	
	/* some default values */
	_debug           = 0;
	verbose          = 1;
	threshold        = 0.001;
	numClusters      = 0;
	isBinaryFile     = 0;
	is_output_timing = 0;
	filename         = NULL;
	center_filename  = NULL;
	output_id        = 1;    
	
	while ( (opt=getopt(argc,argv,"p:i:c:n:f:t:abdohq"))!= EOF) {
		switch (opt) {
			case 'i': filename=optarg;
						break;
			case 'c': center_filename=optarg;
						break;
			case 'b': isBinaryFile = 1;
						break;
			case 't': threshold=atof(optarg);
						break;
			case 'n': numClusters = atoi(optarg);
						break;
			case 'f': output_id = atoi(optarg);
						break;
			case 'o': is_output_timing = 1;
						break;
			case 'q': verbose = 0;
						break;
			case 'd': _debug = 1;
						break;
			default: usage(argv[0], threshold);
						break;
		}
	}
	if (center_filename == NULL)
		center_filename = filename;
	
	if (filename == 0 || numClusters <= 1) usage(argv[0], threshold);
	
	if (is_output_timing) io_timing = wtime();
	
	/* read data points from file ------------------------------------------*/
	printf("reading data points from file %s\n",filename);
	
	objects = file_read(isBinaryFile, filename, &numObjs, &numCoords);
	if (objects == NULL) exit(1);
    
	if (numObjs < numClusters) {
		printf("Error: number of clusters must be larger than the number of data points to be clustered.\n");
		free(objects[0]);
		free(objects);
		return 1;
	}
	//Change the number of coordinates HERE
	numCoords=2;
	/* allocate a 2D space for clusters[] (coordinates of cluster centers)
		this array should be the same across all processes                  */
	clusters    = (double**) malloc(numClusters *             sizeof(double*));
	assert(clusters != NULL);
	clusters[0] = (double*)  malloc(numClusters * numCoords * sizeof(double));
	assert(clusters[0] != NULL);
	for (i=1; i<numClusters; i++)
		clusters[i] = clusters[i-1] + numCoords;
	
	/* read the first numClusters elements from file center_filename as the
		* initial cluster centers*/
	if (center_filename != filename) {
		printf("reading initial %d centers from file %s\n", numClusters,
				center_filename);
		/* read the first numClusters data points from file */
		read_n_objects(isBinaryFile, center_filename, numClusters,
						numCoords, clusters);
	}
	else {
		//Dimension sampling
		//Change dim variable to select the dimension's combination.Select "numCoords" dimensions starting from dimension "dim"
		printf("selecting the first %d elements as initial centers\n",numClusters);
		//Example : numCoords=3,dim=2 => clusters[0][0]=objects[0][2] , clusters[0][1]=objects[0][3] , clusters[0][2]=objects[0][4]
		//                               clusters[1][0]=objects[1][2] , clusters[1][1]=objects[1][3] , clusters[1][2]=objects[1][4]
		//                                                                       ......
		//                               clusters[N][0]=objects[N][2] , clusters[N][1]=objects[N][3] , clusters[N][2]=objects[N][4]
		int dim=3;
		for (j=dim; j<numCoords+dim; j++){
			for (i=0; i<numClusters; i++){
				clusters[i][j-dim] = objects[i][j];
			}
		}
	}
	//Change the number of objects HERE
	numObjs/=1;
	/* check initial cluster centers for repeatition */
	if (check_repeated_clusters(numClusters, numCoords, clusters) == 0) {
		printf("Error: some initial clusters are repeated. Please select distinct initial centers\n");
		return 1;
	}
	
	
	if (is_output_timing) {
		timing            = wtime();
		io_timing         = timing - io_timing;
		clustering_timing = timing;
	}
	
	/* start the timer for the core computation -----------------------------*/
	/* membership: the cluster id for each data object */
	membership = (int*) malloc(numObjs * sizeof(int));
	assert(membership != NULL);
	
	seq_kmeans(objects, numCoords, numObjs, numClusters, threshold, membership,clusters);
	
	free(objects[0]);
	free(objects);
	
	if (is_output_timing) {
		timing            = wtime();
		clustering_timing = timing - clustering_timing;
	}
	
	/* output: the coordinates of the cluster centres ----------------------*/
	fd=fopen("approximated_output.txt","w+");
	if(fd==NULL) {
		printf("Error with fopen().\n");
		exit(1);
	}
	
	if (_debug) {
		printf("Sorted initial cluster centers:\n");
		for (i=0; i<numClusters; i++) {
			printf("clusters[%d]=",i);
			for (j=0; j<numCoords; j++)
				printf(" %6.2f", clusters[i][j]);
			printf("\n");
		}
	}
	file_write("approximated_output.txt", numClusters, numObjs, numCoords, clusters,membership, verbose,output_id);
	
	free(membership);
	free(clusters[0]);
	free(clusters);
	
	/*---- output performance numbers ---------------------------------------*/
	if (is_output_timing) {
		io_timing += wtime() - timing;
		printf(YEL "\nPerforming **** Regular Kmeans (sequential version) ****\n" RESET);
		
		printf("Input file:     %s\n", filename);
		printf("numObjs       = %d\n", numObjs);
		printf("numCoords     = %d\n", numCoords);
		printf("numClusters   = %d\n", numClusters);
		printf("threshold     = %.6f\n", threshold);
		
		printf("I/O time           = %10.4f sec\n", io_timing);
		printf("Computation timing = %10.4f sec\n", clustering_timing);
	}
	
	return(0);
}
