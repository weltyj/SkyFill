# define	NUMBER_OF_CLUSTERS	4
# define	MAXIMUM_ITERATIONS	100
 
 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include "tiffio.h"

// in case Windows MSVC, or 
#define _USE_MATH_DEFINES
#include <math.h>

#include "skyfill_tif.h"
#include "kmeans_rgb.h"
 

typedef struct {
    int	x,y;
    float r,g,b ;
    int		group;
} POINT;


POINT * load_pixels_for_kmeans(tdata_t *image,int num_pts, int *x_list, int *y_list, struct sky_pixel_data *sp_list)
{
    int		i;

    POINT * pts = (POINT*) malloc(sizeof(POINT) * num_pts);

    for ( i = 0; i < num_pts; i++ ) {
	pts[i].x = x_list[i] ;
	pts[i].y = y_list[i] ;
	pts[i].r = sp_list[i].rgb_err[0] ;
	pts[i].g = sp_list[i].rgb_err[1] ;
	pts[i].b = sp_list[i].rgb_err[2] ;
	/*  		pts[i].r = ((uint16_t *)(image[pts[i].y]))[IMAGE_NSAMPLES*pts[i].x+0] ;  */
	/*  		pts[i].g = ((uint16_t *)(image[pts[i].y]))[IMAGE_NSAMPLES*pts[i].x+1] ;  */
	/*  		pts[i].b = ((uint16_t *)(image[pts[i].y]))[IMAGE_NSAMPLES*pts[i].x+2] ;  */
    }
    return pts;	
}

/*-------------------------------------------------------
  dist2

  This function returns the squared euclidean distance
  between two data points.
  -------------------------------------------------------*/
float dist2(POINT * a, POINT * b)
{
    float R = a->r - b->r;
    float G = a->g - b->g;
    float B = a->b - b->b;
    return R*R + G*G + B*B ;
}

/*------------------------------------------------------
  nearest

  This function returns the index of the cluster centroid
  nearest to the data point passed to this function.
  ------------------------------------------------------*/
int nearest(POINT * pt, POINT * cent, int n_cluster)
{
    int i, clusterIndex;
    double d, min_d;

    min_d = HUGE_VAL;
    clusterIndex = pt->group;	
    for (i = 0; i < n_cluster; i++) {
	d = dist2(&cent[i], pt);
	if ( d < min_d ) {
	    min_d = d;
	    clusterIndex = i;
	}
    }	
    return clusterIndex;
}

/*------------------------------------------------------
  nearestDistance

  This function returns the distance of the cluster centroid
  nearest to the data point passed to this function.
  ------------------------------------------------------*/
double nearestDistance(POINT * pt, POINT * cent, int n_cluster)
{
    int i;
    double d, min_d;

    min_d = HUGE_VAL;
    for (i = 0; i < n_cluster; i++) {
	d = dist2(&cent[i], pt);
	if ( d < min_d ) {
	    min_d = d;
	}
    }

    return min_d;
}

/*----------------------------------------------------------------------
  bisectionSearch

  This function makes a bisectional search of an array of values that are
  ordered in increasing order, and returns the index of the first element
  greater than the search value passed as a parameter.

  This code is adapted from code by Andy Allinger given to the public
  domain.

Input:
x	A pointer to an array of values in increasing order to be searched.
n	The number of elements in the input array x.
v	The search value.
Output:
Returns the index of the first element greater than the search value, v.
----------------------------------------------------------------------*/
int bisectionSearch(double *x, int n, double v)
{
    int il, ir, i;


    if (n < 1) {
	return 0;
    }
    /* If v is less than x(0) or greater than x(n-1)  */
    if (v < x[0]) {
	return 0;
    }
    else if (v > x[n-1]) {
	return n - 1;
    }

    /*bisection search */
    il = 0;
    ir = n - 1;

    i = (il + ir) / 2;
    while ( i != il ) {
	if (x[i] <= v) {
	    il = i;
	} else {
	    ir = i;
	}
	i = (il + ir) / 2;		
    }		

    if (x[i] <= v)
	i = ir;
    return i;
} /* end of bisectionSearch */

/*-------------------------------------------------------
  kppAllinger

  This function uses the K-Means++ method to select
  the cluster centroids.

  This code is adapted from code by Andy Allinger given to the
  public domain.

Input:
pts		A pointer to an array of data points.
num_pts		The number of points in the pts array.
centroids	A pointer to an array to receive the centroids.
num_clusters	The number of clusters to be found.

Output:
centroids	A pointer to the array of centroids found.	
-------------------------------------------------------*/
void kppAllinger(POINT * pts, int num_pts, POINT * centroids,
	int num_clusters, int first_cluster_pt_index)
{
    int j;
    int selectedIndex;
    int cluster;
    double sum;
    double d;
    double random;	
    double * cumulativeDistances;
    double * shortestDistance;


    cumulativeDistances = (double*) malloc(sizeof(double) * num_pts);
    shortestDistance = (double*) malloc(sizeof(double) * num_pts);	



    if(first_cluster_pt_index >= 0) {
	// calling routine has specified the point for the values of the first cluster
	selectedIndex = first_cluster_pt_index ;
    } else {
	/* Pick the first cluster centroids at random. */
	selectedIndex = rand() % num_pts;
    }

    /* Pick the first cluster centroids at random. */
    centroids[0] = pts[ selectedIndex ];

    for (j = 0; j < num_pts; ++j)
	shortestDistance[j] = HUGE_VAL;	

    /* Select the centroids for the remaining clusters. */
    for (cluster = 1; cluster < num_clusters; cluster++) {

	/* For each point find its closest distance to any of
	   the previous cluster centers */
	for ( j = 0; j < num_pts; j++ ) {
	    d = dist2(&pts[j], &centroids[cluster-1] );

	    if (d < shortestDistance[j])
		shortestDistance[j] = d;
	}

	/* Create an array of the cumulative distances. */
	sum = 0.0;
	for (j = 0; j < num_pts; j++) {
	    sum += shortestDistance[j];
	    cumulativeDistances[j] = sum;
	}

	/* Select a point at random. Those with greater distances
	   have a greater probability of being selected. */
	random = (float) rand() / (float) RAND_MAX * sum;
	selectedIndex = bisectionSearch(cumulativeDistances, num_pts, random);

	/* assign the selected point as the center */
	centroids[cluster] = pts[selectedIndex];
    }

    /* Assign each point the index of it's nearest cluster centroid. */
    for (j = 0; j < num_pts; j++)
	pts[j].group = nearest(&pts[j], centroids, num_clusters);

    free(shortestDistance);
    free(cumulativeDistances);

    return;
}	/* end, kppAllinger */

/*-------------------------------------------------------
  kpp

  This function uses the K-Means++ method to select
  the cluster centroids.
  -------------------------------------------------------*/
void kpp(POINT * pts, int num_pts, POINT * centroids,
	int num_clusters, int first_cluster_pt_index)
{
    int j;
    int cluster;
    double sum;
    double * distances;


    distances = (double*) malloc(sizeof(double) * num_pts);

    if(first_cluster_pt_index >= 0) {
	// calling routine has specified the point for the values of the first cluster
	centroids[0] = pts[first_cluster_pt_index] ;
    } else {
	/* Pick the first cluster centroids at random. */
	centroids[0] = pts[ rand() % num_pts ];
    }


    /* Select the centroids for the remaining clusters. */
    for (cluster = 1; cluster < num_clusters; cluster++) {

	/* For each data point find the nearest centroid, save its
	   distance in the distance array, then add it to the sum of
	   total distance. */
	sum = 0.0;
	for ( j = 0; j < num_pts; j++ ) {
	    distances[j] = 
		nearestDistance(&pts[j], centroids, cluster);
	    sum += distances[j];
	}

	/* Find a random distance within the span of the total distance. */
	sum = sum * rand() / (RAND_MAX - 1);

	/* Assign the centroids. the point with the largest distance
	   will have a greater probability of being selected. */
	for (j = 0; j < num_pts; j++ ) {
	    sum -= distances[j];
	    if ( sum <= 0)
	    {
		centroids[cluster] = pts[j];
		break;
	    }
	}
    }

    /* Assign each observation the index of it's nearest cluster centroid. */
    for (j = 0; j < num_pts; j++)
	pts[j].group = nearest(&pts[j], centroids, num_clusters);

    free(distances);

    return;
}	/* end, kpp */


/*-------------------------------------------------------
  lloyd

  This function clusters the data using Lloyd's K-Means algorithm
  after selecting the intial centroids using the K-Means++
  method.
  It returns a pointer to the memory it allocates containing
  the array of cluster centroids.
  -------------------------------------------------------*/
POINT * lloyd(POINT * pts, int num_pts, int num_clusters, int maxTimes, int first_cluster_pt_index)
{
    int i, clusterIndex;
    int changes;
    int acceptable = num_pts / 1000;	/* The maximum point changes acceptable. */


    if (num_clusters == 1 || num_pts <= 0 || num_clusters > num_pts )
	return 0;


    POINT * centroids = (POINT *)malloc(sizeof(POINT) * num_clusters);

    if ( maxTimes < 1 )
	maxTimes = 1;

    /*	Assign initial clustering randomly using the Random Partition method
	for (i = 0; i < num_pts; i++ ) {
	pts[i].group = i % num_clusters;
	}
     */

    /* or use the k-Means++ method */

    /* Original version
       kpp(pts, num_pts, centroids, num_clusters, first_cluster_pt_index);
     */
    /* Faster Allinger version */
    kppAllinger(pts, num_pts, centroids, num_clusters, first_cluster_pt_index);

    do {
	/* Calculate the centroid of each cluster.
	   ----------------------------------------*/

	/* Initialize the x, y and cluster totals. */
	for ( i = 0; i < num_clusters; i++ ) {
	    centroids[i].group = 0;		/* used to count the cluster members. */
	    centroids[i].r = 0;			/* used for r value totals. */
	    centroids[i].g = 0;			/* used for g value totals. */
	    centroids[i].b = 0;			/* used for b value totals. */
	}

	/* Add each observation's x and y to its cluster total. */
	for (i = 0; i < num_pts; i++) {
	    clusterIndex = pts[i].group;
	    centroids[clusterIndex].group++;
	    centroids[clusterIndex].r += pts[i].r;
	    centroids[clusterIndex].g += pts[i].g;
	    centroids[clusterIndex].b += pts[i].b;
	}

	/* Divide each cluster's x and y totals by its number of data points. */
	for ( i = 0; i < num_clusters; i++ ) {
	    centroids[i].r /= centroids[i].group;
	    centroids[i].g /= centroids[i].group;
	    centroids[i].b /= centroids[i].group;
	}

	/* Find each data point's nearest centroid */
	changes = 0;
	for ( i = 0; i < num_pts; i++ ) {
	    clusterIndex = nearest(&pts[i], centroids, num_clusters);
	    if (clusterIndex != pts[i].group) {
		pts[i].group = clusterIndex;
		changes++;
	    }
	}

	maxTimes--;
    } while ((changes > acceptable) && (maxTimes > 0));

    /* Set each centroid's group index */
    for ( i = 0; i < num_clusters; i++ )
	centroids[i].group = i;

    return centroids;
}	/* end, lloyd */

void classify_pixels(int n_pixels, POINT * pts, POINT * centroids, int num_clusters, int *classification)
{
    // the cluster with the largest blue value is the sky (hopefully)

    int max_b=0 ;
    int sky_cluster=0 ;
    int sky_group=0 ;
    int sky_group2=0 ;

    for (int i = 0; i < num_clusters; i++) {
	if(centroids[i].b > max_b) {
	    max_b = centroids[i].b ;
	    sky_group = centroids[i].group ;
	    sky_cluster=i ;
	}
    }

    // compute RGB euclidian distance between other clusters and the sky cluster
    // if it is small the other cluster is likely sky pixels too

    int sky_cluster2 = nearest(&centroids[sky_cluster], centroids, num_clusters) ;
    sky_group2=centroids[sky_cluster2].group ;


    for (int i = 0; i < n_pixels; i++) {
	if(pts[i].group == sky_group) {
	    classification[i] |= SKY_PIXEL_TYPE ;
	}
	if(pts[i].group == sky_group2) {
	    classification[i] |= SKY2_PIXEL_TYPE ;
	}
    }
}

void kmeans_analysis(tdata_t *image, int n_pixels, int *x_list, int *y_list, struct sky_pixel_data *sp_list, int *classification, int first_cluster_pt_index)
{
    POINT * pts;
    POINT * centroids;

    int		num_clusters = NUMBER_OF_CLUSTERS;
    int		maxTimes = MAXIMUM_ITERATIONS;

    pts = load_pixels_for_kmeans(image, n_pixels, x_list, y_list, sp_list) ;

    /* Cluster using the Lloyd algorithm and K-Means++ initial centroids. */
    centroids = lloyd(pts, n_pixels, num_clusters, maxTimes, first_cluster_pt_index);
    classify_pixels(n_pixels, pts, centroids, num_clusters, classification) ;
    free(pts) ;
    free(centroids) ;
}
