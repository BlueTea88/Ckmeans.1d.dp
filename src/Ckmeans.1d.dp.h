/*
 Ckmeans_1d_dp.h --- Head file for Ckmeans.1d.dp
 Declare wrap function "kmeans_1d_dp()"

 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Created: Oct 10, 2010
 */

#include <cstddef> // For size_t
#include <vector>
#include <string>

/* One-dimensional cluster algorithm implemented in C++ */
void kmeans_1d_dp(
    const double *x, const size_t N, 
    const double *z, const size_t M,
    const double *w,
    size_t k_in, int* cluster, double* centers, double* size);

void fill_weighted_dp_matrix(
    const std::vector<double> & x,
    const std::vector< std::vector< double > > & z,
    const std::vector< std::vector< double > > & w,
    std::vector< std::vector< double > > & S,
    std::vector< std::vector< size_t > > & J);

void backtrack_weighted(
    const std::vector<double> & x, 
    const std::vector< std::vector< double > > & w,
    const std::vector< std::vector< size_t > > & J,
    int* cluster, double* centers, double* weights);
