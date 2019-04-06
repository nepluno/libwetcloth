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

#ifndef THREAD_UTILS
#define THREAD_UTILS

#include <thread>
#include <tbb/tbb.h>
#include <vector>

// #define DEBUG_PARALLEL 1
// #define NO_PARALLEL 1

namespace threadutils {
	inline unsigned get_num_threads()
	{
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
		return std::thread::hardware_concurrency();
#else
		return 1U;
#endif
	}

	template<typename Index, typename Callable>
	static void for_each(Index start, Index end, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
		tbb::parallel_for(start, end, 1, func);
#else
		for(Index i = start; i < end; ++i) {
			func(i);
		}
#endif
	}
	
	template<typename Data, typename Callable>
	static void for_each(std::vector<Data>& vec, Callable func) {
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
		tbb::parallel_for(0, (int) vec.size(), 1, [&] (int i) {
			func(vec[i]);
		});
#else
		for(auto& data : vec) {
			func(data);
		}
#endif
	}
	
	template<typename Index, typename Data, typename Callable>
	static Data reduction(Data* array, Index n, const Data& base, Callable func)
	{
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
		return tbb::parallel_reduce( tbb::blocked_range<Data*>( array, array + n ), base, [&](const tbb::blocked_range<Data*>& r, Data init) -> Data {
			for( Data* a = r.begin(); a != r.end(); ++a ) {
				init = func(init, *a);
			}
			return init;
		}, func );
#else
		Data res = base;
		for(Index i = 0; i < n; ++i)
		{
			res = func(res, array[i]);
		}
		return res;
#endif
	}
};
#endif /* ThreadUtils_hpp */
