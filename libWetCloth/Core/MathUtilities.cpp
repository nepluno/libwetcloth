//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "MathUtilities.h"

#include <iomanip>

#include "ThreadUtils.h"

namespace mathutils {
bool approxSymmetric(const MatrixXs& A, const scalar& eps) {
  for (int i = 0; i < A.rows(); ++i)
    for (int j = i + 1; j < A.cols(); ++j)
      if (fabs(A(i, j) - A(j, i)) >= eps) return false;
  return true;
}

scalar scalarRand(const scalar min, const scalar max) {
  static thread_local std::mt19937 generator;
  std::uniform_real_distribution<scalar> distribution(min, max);
  return distribution(generator);
}

void print_histogram_analysis(const std::vector<VectorXs>& vec, int num_bins,
                              const std::string& name, bool nonzero) {
  if (!vec.size()) return;

  VectorXs bucket_mins(vec.size());
  bucket_mins.setConstant(std::numeric_limits<scalar>::infinity());

  VectorXs bucket_maxs(vec.size());
  bucket_maxs.setConstant(-std::numeric_limits<scalar>::infinity());

  const int num_buckets = (int)vec.size();

  threadutils::for_each(0, num_buckets, [&](int bucket_idx) {
    if (vec[bucket_idx].size() == 0) return;

    bucket_maxs[bucket_idx] = vec[bucket_idx].maxCoeff();
    bucket_mins[bucket_idx] = vec[bucket_idx].minCoeff();
  });

  scalar total_min = bucket_mins.minCoeff();
  scalar total_max = bucket_maxs.maxCoeff();
  scalar range = (total_max - total_min);

  if (range == 0.0) return;

  std::vector<std::vector<int> > histogram_by_buckets(vec.size());

  threadutils::for_each(0, num_buckets, [&](int bucket_idx) {
    const VectorXs& nodes = vec[bucket_idx];
    std::vector<int>& hist = histogram_by_buckets[bucket_idx];

    hist.resize(num_bins, 0);

    const int num_nodes = (int)nodes.size();
    for (int i = 0; i < num_nodes; ++i) {
      if (nonzero && nodes[i] == 0.0) continue;

      int idx_bin = std::min(num_bins - 1, (int)((nodes[i] - total_min) /
                                                 range * (scalar)num_bins));
      hist[idx_bin]++;
    }
  });

  std::vector<int> total_bins(num_bins, 0);
  threadutils::for_each(0, num_bins, [&](int bin_idx) {
    for (int j = 0; j < num_buckets; ++j) {
      total_bins[bin_idx] += histogram_by_buckets[j][bin_idx];
    }
  });

  int total_samples = 0;
  for (int i = 0; i < num_bins; ++i) total_samples += total_bins[i];

  std::cout << "[Hist. Anal. for " << name << ", min: " << total_min
            << ", max: " << total_max << "]" << std::endl;
  for (int i = 0; i < num_bins; ++i) {
    std::cout << std::setprecision(2)
              << ((scalar)total_bins[i] / (scalar)total_samples * 100.0)
              << "% ";
  }
  std::cout << std::endl;
}
}  // namespace mathutils
