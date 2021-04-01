#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strtok() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     /* read(), close() */
#include <errno.h>
#include <assert.h>
#include <math.h>
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define MAX_NUMBER 10000000.0

#define MAX_CHAR_PER_LINE 128
int _debug;
double** file_read(int   isBinaryFile,  /* flag: 0 or 1 */
					char *filename,      /* input file name */
					int  *numObjs,       /* no. data objects (local) */
					int  *numCoords)     /* no. coordinates */
{
	double **objects;
	int  i,j,len;
	ssize_t numBytesRead;
	
	if (isBinaryFile) {  /* input file is in raw binary format -------------*/
		int infile;
		if ((infile = open(filename, O_RDONLY, "0600")) == -1) {
			fprintf(stderr, "Error: no such file (%s)\n", filename);
			return NULL;
		}
		numBytesRead = read(infile, numObjs,    sizeof(int));
		assert(numBytesRead == sizeof(int));
		numBytesRead = read(infile, numCoords, sizeof(int));
		assert(numBytesRead == sizeof(int));
		if (_debug) {
			printf("File %s numObjs   = %d\n",filename,*numObjs);
			printf("File %s numCoords = %d\n",filename,*numCoords);
		}
		
		/* allocate space for objects[][] and read all objects */
		len = (*numObjs) * (*numCoords);
		objects    = (double**)malloc((*numObjs) * sizeof(double*));
		assert(objects != NULL);
		objects[0] = (double*) malloc(len * sizeof(float));
		assert(objects[0] != NULL);
		for (i=1; i<(*numObjs); i++)
			objects[i] = objects[i-1] + (*numCoords);
		
		numBytesRead = read(infile, objects[0], len*sizeof(float));
		assert(numBytesRead == len*sizeof(float));
		close(infile);
	}
	else{
		FILE *infile;
		char *line, *ret;
		int   lineLen;
		
		if ((infile = fopen(filename, "r")) == NULL) {
			fprintf(stderr, "Error: no such file (%s)\n", filename);
			return NULL;
		}
		
		/* first find the number of objects */
		lineLen = MAX_CHAR_PER_LINE;
		line = (char*) malloc(lineLen);
		assert(line != NULL);
		
		(*numObjs) = 0;
		while (fgets(line, lineLen, infile) != NULL) {
			/* check each line to find the max line length */
			while (strlen(line) == lineLen-1) {
				/* this line read is not complete */
				len = strlen(line);
				fseek(infile, -len, SEEK_CUR);
				
				/* increase lineLen */
				lineLen += MAX_CHAR_PER_LINE;
				line = (char*) realloc(line, lineLen);
				assert(line != NULL);

				ret = fgets(line, lineLen, infile);
				assert(ret != NULL);
			}

			if (strtok(line, " \t\n") != 0)
				(*numObjs)++;
		}
		rewind(infile);

		/* find the no. coordinates of each object */
		(*numCoords) = 0;
		while (fgets(line, lineLen, infile) != NULL) {
			if (strtok(line, " \t\n") != 0) {
				/* ignore the id (first coordiinate): numCoords = 1; */
				while (strtok(NULL, " ,\t\n") != NULL) (*numCoords)++;
				break; /* this makes read from 1st object */
			}
		}
		rewind(infile);
		if (_debug) {
			printf("File %s numObjs   = %d\n",filename,*numObjs);
			printf("File %s numCoords = %d\n",filename,*numCoords);
		}
		
		/* allocate space for objects[][] and read all objects */
		len = (*numObjs) * (*numCoords);
		objects    = (double**)malloc((*numObjs) * sizeof(double*));
		assert(objects != NULL);
		objects[0] = (double*) malloc(len * sizeof(double));
		assert(objects[0] != NULL);
		for (i=1; i<(*numObjs); i++)
			objects[i] = objects[i-1] + (*numCoords);
			
		i = 0;
		/* read all objects */
		while (fgets(line, lineLen, infile) != NULL) {
			if (strtok(line, " \t\n") == NULL) continue;
			for (j=0; j<(*numCoords); j++) {
				objects[i][j] = atof(strtok(NULL, " ,\t\n"));
			}
			i++;
		}
		assert(i == *numObjs);
		
		fclose(infile);
		free(line);
		
	}
	return objects;
}
double euclid_dist_2(int    numdims,  /* no. dimensions */
                    const double *coord1,   /* [numdims] */
                    const double *coord2)   /* [numdims] */
{
	int i;
	double ans=0.0;
	double diff;    
	
	for (i=0; i<numdims; i++){
		diff = coord1[i]-coord2[i]; //Common Subexpression Elimination
		ans += diff * diff;
	}
	ans = sqrt(ans);
	return(ans);
}

double check_accuracy(double **clusters,double ** approximated_clusters,int numClusters,int numCoords,int *cluster_match_array,double *divergence){
	double error_percentage=0.0,min=MAX_NUMBER,min_conflict=MAX_NUMBER;
	int i,k,j,min_position=0,conflict;
	double error=0.0,cluster_i_dist,average_distance=0.0;
	
	//Check all the clusters from the approximated version(approximated_clusters array)
	for(i=0;i<numClusters;i++){
		//Find the minimum euclidean distance between approximated_clusters[i](C'i) and cluster[k](Ck).So,C'i-Ck is min and min_position=k
		for(k=0;k<numClusters;k++){
			cluster_i_dist = euclid_dist_2(numCoords, clusters[k], approximated_clusters[i]);
			if(cluster_i_dist<min){
				min=cluster_i_dist;
				min_position=k;
			}
		}
		//Check if the Ck cluster has been matched with another cluster from the approximated_clusters array.
		for(j=0;j<numClusters;j++){
			if(cluster_match_array[j]==min_position){break;}
		}
		if(j==numClusters){
			//There is no match,continue with the next cluster(i++)
			cluster_match_array[i]=min_position;
			divergence[i]=min;
			min=MAX_NUMBER;
		}
		else{
			//We have conflict.Check if (C'i-Cmin_position)<(C'j-Cmin_position) or if (C'i-Cmin_position)>(C'j-Cmin_position)
			cluster_i_dist = euclid_dist_2(numCoords, clusters[min_position], approximated_clusters[j]);
			if(min<=cluster_i_dist){
				//The distance between cluster "i" and cluster "min_position" is smaller than the distance between the cluster "j(conflict cluster)" and the cluster "min_position"
				//So,the cluster with the conflict is the cluster j and it needs to be redefined.
				cluster_match_array[i]=min_position;
				divergence[i]=min;
				conflict=j;
				min_conflict=min; 
				min=MAX_NUMBER;
			}
			else if(min>cluster_i_dist){
				//The distance between cluster "i" and cluster "min_position" is bigger than the distance between the cluster "j(conflict cluster)" and the cluster "min_position"
				//So,the cluster with the conflict is the cluster i and it needs to be redefined.
				//printf("conflict in cluster %d,cluster_i_dist=%f\n",i,cluster_i_dist);
				min_conflict=min;    
				conflict=i;
				min=MAX_NUMBER;
			}
			while(1){
				//Find the minimum euclidean distance between approximated_clusters[conflict](C'conflict) and cluster[k](Ck).So,C'conflict-Ck is min and min_position=k
				for(k=0;k<numClusters;k++){
					cluster_i_dist = euclid_dist_2(numCoords, clusters[k], approximated_clusters[conflict]);
					if(cluster_i_dist<min && cluster_i_dist>min_conflict){
						min=cluster_i_dist;
						min_position=k;
					}
				}
				if(min==MAX_NUMBER){
					for(k=0;k<numClusters;k++){
						if(divergence[k]==-1){
							divergence[k]=cluster_i_dist = euclid_dist_2(numCoords, clusters[k], approximated_clusters[conflict]);
							cluster_match_array[conflict]=k;
							break;
						}
					}
					break;
				}
				//Check if the Ck cluster has been matched with another cluster from the approximated_clusters array.
				for(j=0;j<numClusters;j++){
					if(cluster_match_array[j]==min_position){break;}
				}
				if(j==numClusters){
					//There is no match,continue with the next cluster(i++)
					cluster_match_array[conflict]=min_position;
					divergence[conflict]=min;
					break;
				}
				else{
					//We have conflict.Check if (C'conflict-Cmin_position)<(C'j-Cmin_position) or if (C'conflict-Cmin_position)>(C'j-Cmin_position)
					cluster_i_dist = euclid_dist_2(numCoords, clusters[min_position], approximated_clusters[j]);
					if(min<=cluster_i_dist){
						//The distance between cluster "conflict" and cluster "min_position" is smaller than the distance between the cluster "j" and the cluster "min_position"
						//So,the cluster with the conflict is the cluster j and it needs to be redefined.
						cluster_match_array[conflict]=min_position;
						divergence[conflict]=min;
						conflict=j;
						min_conflict=cluster_i_dist;   //cluster_i_dist
						min=MAX_NUMBER;
					}
					else if(min>cluster_i_dist){
						//The distance between cluster "conflict" and cluster "min_position" is bigger than the distance between the cluster "j" and the cluster "min_position"
						//So,the cluster with the conflict remains the same and it needs to be redefined.
						min_conflict=min;   
						min=MAX_NUMBER;
					}
				}
			}
			min=MAX_NUMBER;
		}
	}
	for(i=0;i<numClusters;i++){
		
		if(_debug){printf("divergence[%d]=%f\n",i,divergence[i]);}
		if(divergence[i]>0.0001){error+=1.0;}
		average_distance+=divergence[i];
	}
	average_distance/=numClusters;
	printf(ANSI_COLOR_YELLOW"***********************************************************************************************************************\n"ANSI_COLOR_RESET);
	printf("Errors:%6.2f centroids have wrong coordinates\n",error);
	error_percentage=(100*error)/numClusters;
	printf("Average_distance=%f\n",average_distance);
	return(100-error_percentage);
}
int main(int argc,char *argv[]){
	int     opt;
	double accuracy;
	extern char   *optarg;
	extern int     optind;
			int     i;
			int     numCoords=9, numObjs,numClusters;
			double **clusters;      /* [numClusters][numCoords] cluster center */
			double ** approximated_clusters;
			int *cluster_match_array;
			double *divergence;
			
	/* some default values */
	_debug           = 0;
	
	while ( (opt=getopt(argc,argv,"p:i:c:n:f:t:abdohq"))!= EOF) {
		switch (opt) {
			case 'n': numClusters = atoi(optarg);
						break;
			case 'd': _debug = 1;
						break;
			case 'h':
			default: 
						break;
		}
	}
	cluster_match_array=(int*) malloc(numClusters *sizeof(int));
	for (i=0; i<numClusters; i++){
		cluster_match_array[i]=-1;
	}
	divergence=(double*)  malloc(numClusters*sizeof(double));
	for (i=0; i<numClusters; i++){
		divergence[i]=0.0;
	}
	clusters    = (double**) malloc(numClusters *             sizeof(double*));
	assert(clusters != NULL);
	clusters[0] = (double*)  malloc(numClusters * (numCoords) * sizeof(double));
	assert(clusters[0] != NULL);
	for (i=1; i<numClusters; i++){
		clusters[i] = clusters[i-1] + numCoords;
	}
	clusters = file_read(0,"output.txt", &numObjs, &numCoords);
	//Change the number of coordinates HERE
	numCoords=2;
	
	approximated_clusters    = (double**) malloc(numClusters *sizeof(double*));
	assert(approximated_clusters != NULL);
	approximated_clusters[0] = (double*)  malloc(numClusters * numCoords * sizeof(double));
	assert(approximated_clusters[0] != NULL);
	for (i=1; i<numClusters; i++){
		approximated_clusters[i] = approximated_clusters[i-1] + numCoords;
	}
	approximated_clusters = file_read(0,"approximated_output.txt", &numObjs, &numCoords);
	
	accuracy=check_accuracy(clusters,approximated_clusters,numClusters,numCoords,cluster_match_array,divergence);
	printf("accuracy of approximated version is %6.2f%%\n",accuracy);
	printf(ANSI_COLOR_YELLOW"***********************************************************************************************************************\n"ANSI_COLOR_RESET);
	return(0);   
}
