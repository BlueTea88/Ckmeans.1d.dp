/*
 Ckmeans_1d_dp_main.cpp --- wrapper function for "kmeans_1d_dp()"
 Created: Oct 10, 2010
 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu
 Modified:
 March 20, 2014. Joe Song. Removed parameter int *k from the function.
 Added "int *Ks" and "int *nK" to provide a range of the number of clusters
 to search for. Made other changes.
 March 29, 2014. Haizhou Wang. Replaced "int *Ks" and "int *nK" by
 "int *minK" and "int *maxK".
 */

#include "Ckmeans.1d.dp.h"

/*Wrapper function to call kmeans_1d_dp()*/
extern "C" {
  /*
   x_data:   An array containing input data to be clustered.
   x_length: Length of the one dimensional array.
   z_data:   An array containing values used to calculate the variance objective.
   z_length: Number of z values per x observation.
   weight:   Associated cell weights for x and z combinations.
   k_in:     Number of clusters.
   cluster:  An array of cluster IDs for each point in x.
   centers:  An array of centers for each cluster.
   size:     An array of sizes of each cluster.
   */

  void Ckmeans_1d_dp(double *x_data, int* x_length, 
                     double *z_data, int* z_length, double *weight, 
                     int* k_in, int* cluster, 
                     double* centers, double* size)
  {
    // Call C++ version one-dimensional clustering algorithm*/
    kmeans_1d_dp(x_data, (size_t)*x_length, 
                 z_data, (size_t)*z_length, weight, 
                 (size_t)(*k_in), cluster, centers, size);
  }
} // End of extern "C"
