/*
 Ckmeans_1d_dp.cpp -- Performs 1-D k-means by a dynamic programming
 approach that is guaranteed to be optimal.
 Joe Song
 Computer Science Department
 New Mexico State University
 joemsong@cs.nmsu.edu
 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu
 Created: May 19, 2007
 Updated: January 23, 2017
 1. Modified to find optimal clusters based on a separate variable matrix.
 */

#include "Ckmeans.1d.dp.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// #define DEBUG
void kmeans_1d_dp(const double *x, const size_t N, 
                  const double *z, const size_t M,
                  const double *w,
                  size_t k_in, int* cluster, double* centers,
                  double* size)
{
  // NOTE: All vectors in this program is considered starting at position 0.
  
  // Create vectors to process data
  std::vector< std::vector< double > > S( k_in, std::vector<double>(N) );
  std::vector< std::vector< size_t > > J( k_in, std::vector<size_t>(N) ); 
      
  std::vector<double> x_vec(N);
  std::vector< std::vector< double > > z_vec( N, std::vector<double>(M) );
  std::vector< std::vector< double > > w_vec( N, std::vector<double>(M) );
  
  for (size_t i=0; i<N; ++i) x_vec[i] = x[i];
  
  for (size_t i=0; i<N; ++i){
  for (size_t j=0; j<M; ++j){
    z_vec[i][j] = z[i*M + j];
    w_vec[i][j] = w[i*M + j];
  }}
  
  // Fill in dynamic programming matrix 
  fill_weighted_dp_matrix(x_vec, z_vec, w_vec, S, J);
  
  // Backtrack to find the clusters beginning and ending indices
  backtrack_weighted(x_vec, w_vec, J, cluster, centers, size);  
  
} //end of kmeans_1d_dp()
