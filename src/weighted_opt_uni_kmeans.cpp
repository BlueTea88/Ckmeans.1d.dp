// weighted_opt_uni_kmeans.cpp
//
// Joe Song
// Created: May 21, 2016. Extracted from Ckmeans.1d.dp.cpp
// Modified: Dec 4, 2016. Fixed an error in calculating mean_x1
#include "Ckmeans.1d.dp.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

//#define DEBUG

void fill_weighted_dp_matrix(const std::vector<double> & x,
                             const std::vector< std::vector< double > > & z,
                             const std::vector< std::vector< double > > & w,
                             std::vector< std::vector< double > > & S,
                             std::vector< std::vector< size_t > > & J)
  /*
  x: One dimension vector to be clustered, must be sorted.
  z: N x M matrix. z[i][j] is the j-th value for x[i].
  w: N x M matrix of weights.
  S: K x N matrix. S[k][i] is the sum of squares of the distance from each x[i]
  to its cluster mean when there are exactly x[i] is the last point in
  cluster k
  J: K x N backtrack matrix
  NOTE: All vector indices in this program start at position 0.
  */
{
  const int K = S.size();
  const int N = S[0].size();
  const int M = z[0].size();

  for(int k = 0; k < K; ++k) {
    S[k][0] = 0.0;
    J[k][0] = 0;
  }

  // Pre-calculate vectors
  std::vector<double> w_cumm_sum(N, 0.0);  // cummulative sum of weights across x
  std::vector<double> wz2_cumm_sum(N, 0.0);  // cummulative weighted sum of squared z across x
  std::vector< std::vector<double> > w_cumm(N, std::vector<double>(M)); // cummulative w
  std::vector< std::vector<double> > wz_cumm(N, std::vector<double>(M));  // cummulative weighted sum of z

  for(int i = 0; i < N; ++i) {
  for(int j = 0; j < M; ++j) {
    if((j == 0) && (i > 0)) {
      w_cumm_sum[i] = w_cumm_sum[i-1];
      wz2_cumm_sum[i] = wz2_cumm_sum[i-1];
    }
    w_cumm_sum[i] += w[i][j];
    wz2_cumm_sum[i] += w[i][j] * z[i][j] * z[i][j];

    if(i > 0) {
      w_cumm[i][j] += w_cumm[i-1][j] + w[i][j];
      wz_cumm[i][j] += wz_cumm[i-1][j] + w[i][j] * z[i][j];
    } else {
      w_cumm[i][j] = w[i][j];
      wz_cumm[i][j] = w[i][j] * z[i][j];
    }
  }}

  // Initialise for 1 cluster
  double sum_sq_mean_z;  // total weighted squared mean of z

  for(int i = 0; i < N; ++i) {
    sum_sq_mean_z = 0.0;  // initialise sum of squared means

    for (int j = 0; j < M; ++j){
      if (w_cumm[i][j] > 0.0){
        sum_sq_mean_z += (wz_cumm[i][j]) * (wz_cumm[i][j]) / w_cumm[i][j];
      }
    }

    S[0][i] = wz2_cumm_sum[i] - sum_sq_mean_z;
    J[0][0] = 0;
  }

  for(int k = 1; k < K; ++k) {
  for(int i = 1; i < N; ++i) {

    // Skip iterations that dont require calculations
    if(i < k) continue;
    if((k == K-1) && (i < N-1)) continue;  // for largest k, only need to iterate over i == N-1

    // Initialise S and J
    // If optimal preceding clusters, no change to S from adding one obs in its own cluster
    S[k][i] = S[k-1][i-1];
    J[k][i] = i;

    // Loop across alternative break points
    int jlow = k;
    int jhigh = i-1;

    for(int j = jlow; j <= jhigh; ++j) {

      // Calculate the sum of squares from obs j to i
      double ssq_mean = 0.0;  // sum of weighted squared mean z
      double swz = 0.0;  // weighted sum of z from j to i
      double sji = 0.0;  // weighted sum of squares from j to i

      for(int l = 0; l < M; ++l) {
        if((w_cumm[i][l] - w_cumm[j-1][l]) > 0.0) {
          swz = wz_cumm[i][l] - wz_cumm[j-1][l];
          ssq_mean += (swz * swz) / (w_cumm[i][l] - w_cumm[j-1][l]);
        }
      }

      sji = wz2_cumm_sum[i] - wz2_cumm_sum[j-1] - ssq_mean;

      // Check if the alternate break point has a lower sum of squares
      double SSQ_j = S[k-1][j-1] + sji;
      if(SSQ_j < S[k][i]) {
        S[k][i] = SSQ_j;
        J[k][i] = j;
      }
    }
  }}
}


void backtrack_weighted(const std::vector<double> & x,
                        const std::vector< std::vector< double > > & w,
                        const std::vector< std::vector< size_t > > & J,
                        int* cluster, double* centers, double* weights)
{
  const int K = J.size();
  const size_t N = J[0].size();
  const size_t M = w[0].size();
  int k_sum = 0;

  // Sum weights for each x
  std::vector<double> w_sum(N, 0.0);

  for(size_t i = 0; i < N; ++i) {
  for(size_t j = 0; j < M; ++j) {
    w_sum[i] += w[i][j];
  }}

  // Extract information for number of clusters 1 ... K
  for(int nk = 1; nk <= K; ++nk){

    size_t cluster_right = N-1;
    size_t cluster_left;
    k_sum += nk - 1;  // record the sum of number of clusters so far

    // Backtrack the clusters from the dynamic programming matrix
    for(int k = nk-1; k >= 0; --k) {
      cluster_left = J[k][cluster_right];

      // Record cluster for the ith data point when there are nk clusters
      for(size_t i = cluster_left; i <= cluster_right; ++i)
        cluster[(nk-1)*N + i] = k + 1;  // output clusters as 1-based (not 0-based)

      double sum = 0.0;
      double weight = 0.0;

      for(size_t i = cluster_left; i <= cluster_right; ++i) {
        sum += x[i] * w_sum[i];
        weight += w_sum[i];
      }

      centers[k_sum + k] = sum / weight;
      weights[k_sum + k] = weight;

      if(k > 0) cluster_right = cluster_left - 1;
    }
  }
}
