//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "MathUtilities.h"
#include "ThreadUtils.h"

#include <iomanip>

namespace mathutils
{
    bool approxSymmetric( const MatrixXs& A, const scalar& eps )
    {
        for( int i = 0; i < A.rows(); ++i ) for( int j = i+1; j < A.cols(); ++j ) if( fabs(A(i,j)-A(j,i)) >= eps ) return false;
        return true;
    }
    
    scalar scalarRand( const scalar min, const scalar max )
    {
        static thread_local std::mt19937 generator;
        std::uniform_real_distribution<scalar> distribution(min, max);
        return distribution(generator);
    }
	
	void print_histogram_analysis(const std::vector< VectorXs >& vec, int num_bins, const std::string& name, bool nonzero)
	{
		if(!vec.size()) return;
		
		VectorXs bucket_mins(vec.size());
		bucket_mins.setConstant(std::numeric_limits<scalar>::infinity());
		
		VectorXs bucket_maxs(vec.size());
		bucket_maxs.setConstant(-std::numeric_limits<scalar>::infinity());
		
		const int num_buckets = (int) vec.size();
		
		threadutils::for_each(0, num_buckets, [&] (int bucket_idx) {
			if(vec[bucket_idx].size() == 0) return;
			
			bucket_maxs[bucket_idx] = vec[bucket_idx].maxCoeff();
			bucket_mins[bucket_idx] = vec[bucket_idx].minCoeff();
		});
		
		scalar total_min = bucket_mins.minCoeff();
		scalar total_max = bucket_maxs.maxCoeff();
		scalar range = (total_max - total_min);
		
		if(range == 0.0) return;
		
		std::vector< std::vector< int > > histogram_by_buckets(vec.size());
		
		threadutils::for_each(0, num_buckets, [&] (int bucket_idx) {
			const VectorXs& nodes = vec[bucket_idx];
			std::vector<int>& hist = histogram_by_buckets[bucket_idx];
			
			hist.resize(num_bins, 0);
			
			const int num_nodes = (int) nodes.size();
			for(int i = 0; i < num_nodes; ++i) {
				if(nonzero && nodes[i] == 0.0) continue;
				
				int idx_bin = std::min(num_bins - 1, (int) ((nodes[i] - total_min) / range * (scalar) num_bins));
				hist[idx_bin]++;
			}
		});
		
		std::vector<int> total_bins(num_bins, 0);
		threadutils::for_each(0, num_bins, [&] (int bin_idx) {
			for(int j = 0; j < num_buckets; ++j) {
				total_bins[bin_idx] += histogram_by_buckets[j][bin_idx];
			}
		});
		
		int total_samples = 0;
		for(int i = 0; i < num_bins; ++i) total_samples += total_bins[i];
		
		std::cout << "[Hist. Anal. for " << name <<", min: " << total_min << ", max: " << total_max << "]" << std::endl;
		for(int i = 0; i < num_bins; ++i) {
			std::cout << std::setprecision(2) << ((scalar) total_bins[i] / (scalar) total_samples * 100.0) << "% ";
		}
		std::cout << std::endl;
	}
}

