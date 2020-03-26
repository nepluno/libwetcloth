//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef THREAD_UTILS_H
#define THREAD_UTILS_H

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
	for (Index i = start; i < end; ++i) {
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
	for (auto& data : vec) {
		func(data);
	}
#endif
}

template<typename Index, typename Data, typename Callable>
static Data reduction(Data* array, Index n, const Data& base, Callable func)
{
#if (defined(NDEBUG) || DEBUG_PARALLEL) && !NO_PARALLEL
	return tbb::parallel_reduce( tbb::blocked_range<Data*>( array, array + n ), base, [&](const tbb::blocked_range<Data*>& r, Data init) -> Data {
		for ( Data* a = r.begin(); a != r.end(); ++a ) {
			init = func(init, *a);
		}
		return init;
	}, func );
#else
	Data res = base;
	for (Index i = 0; i < n; ++i)
	{
		res = func(res, array[i]);
	}
	return res;
#endif
}
};
#endif /* ThreadUtils_hpp */
