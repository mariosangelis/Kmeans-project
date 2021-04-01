/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_kmeans.c  (sequential version)                        */
/*   Description:  Implementation of simple k-means clustering algorithm     */
/*                 This program takes an array of N data objects, each with  */
/*                 M coordinates and performs a k-means clustering given a   */
/*                 user-provided value of the number of clusters (K). The    */
/*                 clustering results are saved in 2 arrays:                 */
/*                 1. a returned array of size [K][N] indicating the center  */
/*                    coordinates of K clusters                              */
/*                 2. membership[N] stores the cluster center ids, each      */
/*                    corresponding to the cluster a data object is assigned */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department, Northwestern University                        */
/*            email: wkliao@ece.northwestern.edu                             */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>

#include "kmeans.h"

/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__inline static
float euclid_dist_2(int    numdims,  /* no. dimensions */
                    float *coord1,   /* [numdims] */
                    float *coord2)   /* [numdims] */
{
	int i;
	register float ans=0.0;
	register float diff;    
	
	for (i=0; i<numdims; i++)
		diff = coord1[i]-coord2[i]; //Common Subexpression Elimination
		ans += diff * diff;
	
	return(ans);
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
							int     numCoords,   /* no. coordinates */
							float  *object,      /* [numCoords] */
							float **clusters)    /* [numClusters][numCoords] */
{
	int   index, i;
	float dist, min_dist;
	
	/* find the cluster id that has min distance to object */
	index    = 0;
	min_dist = euclid_dist_2(numCoords, object, clusters[0]);
	
	
	for (i=1; i<numClusters; i++) {
		dist = euclid_dist_2(numCoords, object, clusters[i]);
		/* no need square root */
		if (dist < min_dist) { /* find the min and its array index */
			min_dist = dist;
			index    = i;
		}
	}
	return(index);
}

/*----< seq_kmeans() >-------------------------------------------------------*/
/* return an array of cluster centers of size [numClusters][numCoords]       */
int seq_kmeans(float **objects,      /* in: [numObjs][numCoords] */
				int     numCoords,    /* no. features */
				int     numObjs,      /* no. objects */
				int     numClusters,  /* no. clusters */
				float   threshold,    /* % objects change membership */
				int    *membership,   /* out: [numObjs] */
				float **clusters)     /* out: [numClusters][numCoords] */
				
{
	int      i, j, index, loop=0;
	int     *newClusterSize; /* [numClusters]: no. objects assigned in each
								new cluster */
	float    delta;          /* % of objects change their clusters */
	float  **newClusters;   /* [numClusters][numCoords] */
	int num_of_reps = 0;    //number of reps in order for the system to converge
	int num_of_changes=0;   //number of points that change cluster for every round
	/* initialize membership[] */
	register int temp=0;
	for (i=0; i<numObjs; i++) membership[i] = -1;   //initialize all object to the default cluster with id -1

	/* need to initialize newClusterSize and newClusters[0] to all 0 */
	newClusterSize = (int*) calloc(numClusters, sizeof(int));
	assert(newClusterSize != NULL);

	newClusters    = (float**) malloc(numClusters *            sizeof(float*));
	assert(newClusters != NULL);
	newClusters[0] = (float*)  calloc(numClusters * numCoords, sizeof(float));
	assert(newClusters[0] != NULL);
	for (i=1; i<numClusters; i++)
		newClusters[i] = newClusters[i-1] + numCoords;
		

	do {
		num_of_changes=0;
		num_of_reps++;
		delta = 0.0;
		for (i=0; i<numObjs; i++) {
			/* find the array index of nestest cluster center */
			index = find_nearest_cluster(numClusters, numCoords, objects[i],
											clusters);
			
			//printf("Index:%d for object(%d)\n",index,i);
			/* if membership changes, increase delta by 1 */
			if (membership[i] != index){
				delta += 1.0;       //if the previus cluster id for this point was different change the delta value by 1
				num_of_changes+=1;       
			}
			/* assign the membership to object i */
			membership[i] = index;                          //for the point i, assign to cluster #index 

			/* update new cluster center : sum of objects located within */
			newClusterSize[index]++;                        //the number of points for cluster #index is increased by 1
			for (j=0; j<numCoords; j++)
				newClusters[index][j] += objects[i][j];     //add for every  of the cluster #index 
		}
		//printf("\n");
		/* average the sum and replace old cluster center with newClusters */
		for (i=0; i<numClusters; i++) {
			temp = newClusterSize[i];   //Loop Invariant code motion
			for (j=0; j<numCoords; j++) {
				if (temp > 0)
					clusters[i][j] = newClusters[i][j] / temp;
				newClusters[i][j] = 0.0;   /* set back to 0 */
			}
			newClusterSize[i] = 0;   /* set back to 0 */ 
		}
		
		delta /= numObjs;
	} while (delta > threshold && loop++ < 500);

	printf(YEL "*********************************************\n" RESET);
	printf("Need %d repetitions in order to converge.\n",num_of_reps);
	printf("%d points change clusters in the last round.\n",num_of_changes);
	printf(YEL "*********************************************\n" RESET);

	free(newClusters[0]);
	free(newClusters);
	free(newClusterSize);

	return 1;
}

