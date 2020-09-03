//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <vector>

#include "MathDefs.h"
#include "MathUtilities.h"
#include "ThreadUtils.h"

#ifndef SORTER_H
#define SORTER_H

class Sorter {
 public:
  Sorter();
  Sorter(int ni_, int nj_, int nk_);
  void resize(int ni_, int nj_, int nk_);
  ~Sorter();

  inline bool has_bucket(int i, int j, int k) const {
    return !(i < 0 || i > ni - 1 || j < 0 || j > nj - 1 || k < 0 || k > nk - 1);
  }

  inline bool has_bucket(const Vector3i& h) const {
    return !(h(0) < 0 || h(0) > ni - 1 || h(1) < 0 || h(1) > nj - 1 ||
             h(2) < 0 || h(2) > nk - 1);
  }

  inline bool has_bucket_x(int i, int j, int k) const {
    return !(i < 0 || i > ni || j < 0 || j > nj - 1 || k < 0 || k > nk - 1);
  }

  inline bool has_bucket_x(const Vector3i& h) const {
    return !(h(0) < 0 || h(0) > ni || h(1) < 0 || h(1) > nj - 1 || h(2) < 0 ||
             h(2) > nk - 1);
  }

  inline bool has_bucket_y(int i, int j, int k) const {
    return !(i < 0 || i > ni - 1 || j < 0 || j > nj || k < 0 || k > nk - 1);
  }

  inline bool has_bucket_y(const Vector3i& h) const {
    return !(h(0) < 0 || h(0) > ni - 1 || h(1) < 0 || h(1) > nj || h(2) < 0 ||
             h(2) > nk - 1);
  }

  inline bool has_bucket_z(int i, int j, int k) const {
    return !(i < 0 || i > ni - 1 || j < 0 || j > nj - 1 || k < 0 || k > nk);
  }

  inline bool has_bucket_z(const Vector3i& h) const {
    return !(h(0) < 0 || h(0) > ni - 1 || h(1) < 0 || h(1) > nj - 1 ||
             h(2) < 0 || h(2) > nk);
  }

  inline int bucket_index(int i, int j, int k) const {
    return k * (ni * nj) + j * ni + i;
  }

  inline int bucket_index(const Vector3i& h) const {
    return h(2) * (ni * nj) + h(1) * ni + h(0);
  }

  inline uint64_t hash(int i, int j, int k, int pidx) const {
    int high_part = bucket_index(i, j, k);
    return (uint64_t)high_part << 32UL | (uint64_t)pidx;
  }

  inline Vector3i bucket_handle(int bucket_idx) const {
    Vector3i p;
    p[2] = bucket_idx / (ni * nj);
    p[1] = (bucket_idx - p[2] * (ni * nj)) / ni;
    p[0] = bucket_idx - (p[2] * (ni * nj) + p[1] * ni);

    return p;
  }

  inline Vector3i bucket_handle_x(int bucket_idx) const {
    Vector3i p;
    p[2] = bucket_idx / ((ni + 1) * nj);
    p[1] = (bucket_idx - p[2] * ((ni + 1) * nj)) / (ni + 1);
    p[0] = bucket_idx - (p[2] * ((ni + 1) * nj) + p[1] * (ni + 1));

    return p;
  }

  inline Vector3i bucket_handle_y(int bucket_idx) const {
    Vector3i p;
    p[2] = bucket_idx / (ni * (nj + 1));
    p[1] = (bucket_idx - p[2] * (ni * (nj + 1))) / ni;
    p[0] = bucket_idx - (p[2] * (ni * (nj + 1)) + p[1] * ni);

    return p;
  }

  inline Vector3i bucket_handle_z(int bucket_idx) const {
    Vector3i p;
    p[2] = bucket_idx / (ni * nj);
    p[1] = (bucket_idx - p[2] * (ni * nj)) / ni;
    p[0] = bucket_idx - (p[2] * (ni * nj) + p[1] * ni);

    return p;
  }

  inline Vector3i bucket_handle_with_corner(int bucket_idx) const {
    Vector3i p;
    p[2] = bucket_idx / ((ni + 1) * (nj + 1));
    p[1] = (bucket_idx - p[2] * ((ni + 1) * (nj + 1))) / (ni + 1);
    p[0] = bucket_idx - (p[2] * ((ni + 1) * (nj + 1)) + p[1] * (ni + 1));

    return p;
  }

  inline int dim_size(int d) const { return Vector3i(ni, nj, nk)(d); }

  inline int size() const { return ni * nj * nk; }

  template <typename Callable>
  void sort(size_t total_size, Callable func) {
    if (array_idx.size() != total_size) {
      array_idx.resize(total_size);
    }

    memset(&array_sup[0], 0, array_sup.size() * sizeof(std::pair<int, int>));

    const int np = (int)total_size;

    threadutils::for_each(0, np, [&](int pidx) {
      int i, j, k;
      func(pidx, i, j, k);
      i = std::max(0, std::min(ni - 1, i));
      j = std::max(0, std::min(nj - 1, j));
      k = std::max(0, std::min(nk - 1, k));
      array_idx[pidx] = hash(i, j, k, pidx);
    });

    tbb::parallel_sort(array_idx.begin(), array_idx.end());

    threadutils::for_each(0, np, [&](int pidx) {
      int G_ID = pidx;
      int G_ID_PREV = G_ID - 1;
      int G_ID_NEXT = G_ID + 1;

      unsigned int cell = (unsigned int)(array_idx[G_ID] >> 32UL);
      unsigned int cell_prev =
          G_ID_PREV < 0 ? -1U : (unsigned int)(array_idx[G_ID_PREV] >> 32UL);
      unsigned int cell_next =
          G_ID_NEXT >= np ? -1U : (unsigned int)(array_idx[G_ID_NEXT] >> 32UL);
      if (cell != cell_prev) {
        // I'm the start of a cell
        array_sup[cell].first = G_ID;
      }
      if (cell != cell_next) {
        // I'm the end of a cell
        array_sup[cell].second = G_ID + 1;
      }
    });
  }

  inline int get_bucket_size(int bucket_idx) const {
    if (bucket_idx < 0 || bucket_idx >= (int)array_sup.size()) return -1;

    return array_sup[bucket_idx].second - array_sup[bucket_idx].first;
  }

  template <typename Callable>
  void get_buckets(int i, int j, int k, int wl, int wh, int hl, int hh, int dl,
                   int dh, Callable func) const {
    for (int sk = k + dl; sk <= k + dh; ++sk)
      for (int si = i + wl; si <= i + wh; si++)
        for (int sj = j + hl; sj <= j + hh; sj++) {
          get_bucket(si, sj, sk, func);
        }
  }

  template <typename Callable>
  void fast_sweep_buckets(int dir, Callable func) const {
    const int swCount = dir + 1;
    const int totalLevels = ni + nj + nk;
    const int start =
        (swCount == 2 || swCount == 5 || swCount == 7 || swCount == 8)
            ? totalLevels
            : 3;
    const int end = (start == 3) ? (totalLevels + 1) : 2;
    const int incr = (start == 3) ? true : false;

    const int xSweepOff = (swCount == 4 || swCount == 8) ? ni + 1 : 0;
    const int ySweepOff = (swCount == 2 || swCount == 6) ? nj + 1 : 0;
    const int zSweepOff = (swCount == 3 || swCount == 7) ? nk + 1 : 0;

    for (int level = start; level != end;
         level = (incr) ? level + 1 : level - 1) {
      const int xs = std::max(1, level - (nj + nk));
      const int ys = std::max(1, level - (ni + nk));
      const int xe = std::min(ni, level - 2);
      const int ye = std::min(nj, level - 2);

      const int xr = xe - xs + 1;
      const int yr = ye - ys + 1;

      int tth = xr * yr;  // Total number of threads needed

      threadutils::for_each(0, tth, [&](int tidx) {
        const int bidy = tidx / xr;
        const int bidx = tidx - bidy * xr;

        const int x = bidx + xs;
        const int y = bidy + ys;

        if (x <= ni && y <= nj) {
          const int z = level - (x + y);

          if (z > 0 && z <= nk) {
            const int k = abs(z - zSweepOff);
            const int j = abs(y - ySweepOff);
            const int i = abs(x - xSweepOff);

            const int bucket_idx = bucket_index(i - 1, j - 1, k - 1);
            func(bucket_idx);
          }
        }
      });
    }
  }

  bool empty(int bucket_idx) const {
    const std::pair<int, int>& G_START_END = array_sup[bucket_idx];
    return G_START_END.second == G_START_END.first;
  }

  template <typename Callable>
  void get_bucket(int bucket_idx, Callable func) const {
    const std::pair<int, int>& G_START_END = array_sup[bucket_idx];
    for (int N_ID = G_START_END.first; N_ID < G_START_END.second; ++N_ID) {
      func((int)(array_idx[N_ID] & 0xFFFFFFFFUL));
    }
  }

  template <typename Callable>
  void get_bucket(int i, int j, int k, Callable func) const {
    if (i < 0 || i > ni - 1 || j < 0 || j > nj - 1 || k < 0 || k > nk - 1)
      return;
    int bucket_idx = bucket_index(i, j, k);
    get_bucket(bucket_idx, func);
  }

  template <typename Callable>
  void for_each_bucket(Callable func, bool forced_serial = false) const {
    int nsystem = ni * nj * nk;
    if (forced_serial) {
      for (int bucket_idx = 0; bucket_idx < nsystem; ++bucket_idx) {
        func(bucket_idx);
      }
    } else {
      threadutils::for_each(0, nsystem,
                            [&](int bucket_idx) { func(bucket_idx); });
    }
  }

  template <typename Callable>
  void for_each_bucket_x(Callable func, bool forced_serial = false) const {
    int nsystem = (ni + 1) * nj * nk;
    if (forced_serial) {
      for (int bucket_idx = 0; bucket_idx < nsystem; ++bucket_idx) {
        func(bucket_idx);
      }
    } else {
      threadutils::for_each(0, nsystem,
                            [&](int bucket_idx) { func(bucket_idx); });
    }
  }

  template <typename Callable>
  void for_each_bucket_y(Callable func, bool forced_serial = false) const {
    int nsystem = ni * (nj + 1) * nk;
    if (forced_serial) {
      for (int bucket_idx = 0; bucket_idx < nsystem; ++bucket_idx) {
        func(bucket_idx);
      }
    } else {
      threadutils::for_each(0, nsystem,
                            [&](int bucket_idx) { func(bucket_idx); });
    }
  }

  template <typename Callable>
  void for_each_bucket_z(Callable func, bool forced_serial = false) const {
    int nsystem = ni * nj * (nk + 1);
    if (forced_serial) {
      for (int bucket_idx = 0; bucket_idx < nsystem; ++bucket_idx) {
        func(bucket_idx);
      }
    } else {
      threadutils::for_each(0, nsystem,
                            [&](int bucket_idx) { func(bucket_idx); });
    }
  }

  template <typename Callable>
  void for_each_bucket_with_corner(Callable func,
                                   bool forced_serial = false) const {
    int nsystem = (ni + 1) * (nj + 1) * (nk + 1);
    if (forced_serial) {
      for (int bucket_idx = 0; bucket_idx < nsystem; ++bucket_idx) {
        func(bucket_idx);
      }
    } else {
      threadutils::for_each(0, nsystem,
                            [&](int bucket_idx) { func(bucket_idx); });
    }
  }

  template <typename Callable>
  void for_each_bucket_colored(Callable func, int numcolors = 2) const {
    const int sni = (ni + numcolors - 1) / numcolors;
    const int snj = (nj + numcolors - 1) / numcolors;
    const int snk = (nk + numcolors - 1) / numcolors;

    const int nsystem = sni * snj * snk;

    for (int t = 0; t < numcolors; ++t)
      for (int s = 0; s < numcolors; ++s)
        for (int r = 0; r < numcolors; ++r) {
          threadutils::for_each(0, nsystem, [&](int scaled_bucket_idx) {
            const int sk = scaled_bucket_idx / (sni * snj);
            const int sj = (scaled_bucket_idx - sk * sni * snj) / sni;
            const int si = scaled_bucket_idx - sk * sni * snj - sj * sni;

            const int k = sk * numcolors + t;
            if (k >= nk) return;

            const int j = sj * numcolors + s;
            if (j >= nj) return;

            const int i = si * numcolors + r;
            if (i >= ni) return;

            const int bucket_idx = k * ni * nj + j * ni + i;
            func(bucket_idx);
          });
        }
  }

  template <typename Callable>
  void loop_neighbor_bucket_particles(int bucket_idx, Callable func,
                                      int search_range = 1) const {
    Vector3i center_handle = bucket_handle(bucket_idx);
    for (int i = center_handle(0) - search_range;
         i <= center_handle(0) + search_range; ++i)
      for (int j = center_handle(1) - search_range;
           j <= center_handle(1) + search_range; ++j)
        for (int k = center_handle(2) - search_range;
             k <= center_handle(2) + search_range; ++k) {
          if (i < 0 || i > ni - 1 || j < 0 || j > nj - 1 || k < 0 || k > nk - 1)
            continue;

          int bucket_idx = bucket_index(i, j, k);
          const std::pair<int, int>& G_START_END = array_sup[bucket_idx];
          for (int N_ID = G_START_END.first; N_ID < G_START_END.second;
               ++N_ID) {
            if (func((int)(array_idx[N_ID] & 0xFFFFFFFFUL), bucket_idx)) return;
          }
        }
  }

  template <typename Callable>
  void loop_neighbor_buckets(int bucket_idx, Callable func,
                             int search_range = 1) const {
    Vector3i center_handle = bucket_handle(bucket_idx);
    for (int i = center_handle(0) - search_range;
         i <= center_handle(0) + search_range; ++i)
      for (int j = center_handle(1) - search_range;
           j <= center_handle(1) + search_range; ++j)
        for (int k = center_handle(2) - search_range;
             k <= center_handle(2) + search_range; ++k) {
          if (i < 0 || i > ni - 1 || j < 0 || j > nj - 1 || k < 0 || k > nk - 1)
            continue;

          int bucket_idx = bucket_index(i, j, k);
          if (func(bucket_idx)) return;
        }
  }

  template <typename Callable>
  void for_each_bucket_particles(Callable func) const {
    int nsystem = ni * nj * nk;
    threadutils::for_each(0, nsystem, [&](int bucket_idx) {
      const std::pair<int, int>& G_START_END = array_sup[bucket_idx];
      for (int N_ID = G_START_END.first; N_ID < G_START_END.second; ++N_ID) {
        func((int)(array_idx[N_ID] & 0xFFFFFFFFUL), bucket_idx);
      }
    });
  }

  template <typename Callable>
  void for_each_bucket_particles_colored(Callable func,
                                         int numcolors = 2) const {
    const int sni = (ni + numcolors - 1) / numcolors;
    const int snj = (nj + numcolors - 1) / numcolors;
    const int snk = (nk + numcolors - 1) / numcolors;

    const int nsystem = sni * snj * snk;

    for (int t = 0; t < numcolors; ++t)
      for (int s = 0; s < numcolors; ++s)
        for (int r = 0; r < numcolors; ++r) {
          threadutils::for_each(0, nsystem, [&](int scaled_bucket_idx) {
            const int sk = scaled_bucket_idx / (sni * snj);
            const int sj = (scaled_bucket_idx - sk * sni * snj) / sni;
            const int si = scaled_bucket_idx - sk * sni * snj - sj * sni;

            const int k = sk * numcolors + t;
            if (k >= nk) return;

            const int j = sj * numcolors + s;
            if (j >= nj) return;

            const int i = si * numcolors + r;
            if (i >= ni) return;

            const int bucket_idx = k * ni * nj + j * ni + i;
            const std::pair<int, int>& G_START_END = array_sup[bucket_idx];
            for (int N_ID = G_START_END.first; N_ID < G_START_END.second;
                 ++N_ID) {
              func((int)(array_idx[N_ID] & 0xFFFFFFFFUL), bucket_idx);
            }
          });
        }
  }

  template <typename Callable>
  void for_each_bucket_particles_colored_randomized(Callable func,
                                                    int numcolors = 2) const {
    const int sni = (ni + numcolors - 1) / numcolors;
    const int snj = (nj + numcolors - 1) / numcolors;
    const int snk = (nk + numcolors - 1) / numcolors;

    const int nsystem = sni * snj * snk;

    std::vector<int> rand_vec_t;
    std::vector<int> rand_vec_s;
    std::vector<int> rand_vec_r;

    mathutils::fisherYates(numcolors, rand_vec_t);
    mathutils::fisherYates(numcolors, rand_vec_s);
    mathutils::fisherYates(numcolors, rand_vec_r);

    for (int t : rand_vec_t)
      for (int s : rand_vec_s)
        for (int r : rand_vec_r) {
          threadutils::for_each(0, nsystem, [&](int scaled_bucket_idx) {
            const int sk = scaled_bucket_idx / (sni * snj);
            const int sj = (scaled_bucket_idx - sk * sni * snj) / sni;
            const int si = scaled_bucket_idx - sk * sni * snj - sj * sni;

            const int k = sk * numcolors + t;
            if (k >= nk) return;

            const int j = sj * numcolors + s;
            if (j >= nj) return;

            const int i = si * numcolors + r;
            if (i >= ni) return;

            const int bucket_idx = k * ni * nj + j * ni + i;
            const std::pair<int, int>& G_START_END = array_sup[bucket_idx];
            for (int N_ID = G_START_END.first; N_ID < G_START_END.second;
                 ++N_ID) {
              func((int)(array_idx[N_ID] & 0xFFFFFFFFUL), bucket_idx);
            }
          });
        }
  }

  std::vector<uint64_t> array_idx;
  std::vector<std::pair<int, int> > array_sup;

  int ni;
  int nj;
  int nk;
};

#endif
