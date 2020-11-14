#include "MakeLevelSet3.h"

#include "MathUtilities.h"

//
// This file is part of the libWetCloth open source project
// The MIT License (MIT)
//
// Copyright (c) 2015, Christopher Batty
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// find distance x0 is from segment x1-x2
static scalar point_segment_distance(const Vector3s &x0, const Vector3s &x1,
                                     const Vector3s &x2) {
  Vector3s dx(x2 - x1);
  scalar m2 = dx.squaredNorm();
  // find parameter value of closest point on segment
  scalar s12 = (x2 - x0).dot(dx) / m2;
  if (s12 < 0) {
    s12 = 0;
  } else if (s12 > 1) {
    s12 = 1;
  }
  // and find the distance
  return (x0 - (s12 * x1 + (1 - s12) * x2)).norm();
}

// find distance x0 is from triangle x1-x2-x3
static scalar point_triangle_distance(const Vector3s &x0, const Vector3s &x1,
                                      const Vector3s &x2, const Vector3s &x3) {
  // first find barycentric coordinates of closest point on infinite plane
  Vector3s x13(x1 - x3), x23(x2 - x3), x03(x0 - x3);
  scalar m13 = x13.squaredNorm(), m23 = x23.squaredNorm(), d = x13.dot(x23);
  scalar invdet = 1. / std::max(m13 * m23 - d * d, 1e-30);
  scalar a = x13.dot(x03), b = x23.dot(x03);
  // the barycentric coordinates themselves
  scalar w23 = invdet * (m23 * a - d * b);
  scalar w31 = invdet * (m13 * b - d * a);
  scalar w12 = 1 - w23 - w31;
  if (w23 >= 0 && w31 >= 0 && w12 >= 0) {  // if we're inside the triangle
    return (x0 - (w23 * x1 + w31 * x2 + w12 * x3)).norm();
  } else {        // we have to clamp to one of the edges
    if (w23 > 0)  // this rules out edge 2-3 for us
      return std::min(point_segment_distance(x0, x1, x2),
                      point_segment_distance(x0, x1, x3));
    else if (w31 > 0)  // this rules out edge 1-3
      return std::min(point_segment_distance(x0, x1, x2),
                      point_segment_distance(x0, x2, x3));
    else  // w12 must be >0, ruling out edge 1-2
      return std::min(point_segment_distance(x0, x1, x3),
                      point_segment_distance(x0, x2, x3));
  }
}

static void check_neighbour(const std::vector<Vector3i> &tri,
                            const std::vector<Vector3s> &x, Array3d &phi,
                            Array3i &closest_tri, const Vector3s &gx, int i0,
                            int j0, int k0, int i1, int j1, int k1) {
  if (closest_tri(i1, j1, k1) >= 0) {
    const Vector3i &pqr = tri[closest_tri(i1, j1, k1)];
    scalar d = point_triangle_distance(gx, x[pqr(0)], x[pqr(1)], x[pqr(2)]);
    if (d < phi(i0, j0, k0)) {
      phi(i0, j0, k0) = d;
      closest_tri(i0, j0, k0) = closest_tri(i1, j1, k1);
    }
  }
}

static void sweep(const std::vector<Vector3i> &tri,
                  const std::vector<Vector3s> &x, Array3d &phi,
                  Array3i &closest_tri, const Vector3s &origin, scalar dx,
                  int di, int dj, int dk) {
  int i0, i1;
  if (di > 0) {
    i0 = 1;
    i1 = phi.ni;
  } else {
    i0 = phi.ni - 2;
    i1 = -1;
  }
  int j0, j1;
  if (dj > 0) {
    j0 = 1;
    j1 = phi.nj;
  } else {
    j0 = phi.nj - 2;
    j1 = -1;
  }
  int k0, k1;
  if (dk > 0) {
    k0 = 1;
    k1 = phi.nk;
  } else {
    k0 = phi.nk - 2;
    k1 = -1;
  }
  for (int k = k0; k != k1; k += dk)
    for (int j = j0; j != j1; j += dj)
      for (int i = i0; i != i1; i += di) {
        Vector3s gx(i * dx + origin[0], j * dx + origin[1], k * dx + origin[2]);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j, k);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i, j - dj, k);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j - dj,
                        k);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i, j, k - dk);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j,
                        k - dk);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i, j - dj,
                        k - dk);
        check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j - dj,
                        k - dk);
      }
}

// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate
// triangle)
static int orientation(scalar x1, scalar y1, scalar x2, scalar y2,
                       scalar &twice_signed_area) {
  twice_signed_area = y1 * x2 - x1 * y2;
  if (twice_signed_area > 0)
    return 1;
  else if (twice_signed_area < 0)
    return -1;
  else if (y2 > y1)
    return 1;
  else if (y2 < y1)
    return -1;
  else if (x1 > x2)
    return 1;
  else if (x1 < x2)
    return -1;
  else
    return 0;  // only true when x1==x2 and y1==y2
}

// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
// if true is returned, the barycentric coordinates are set in a,b,c.
static bool point_in_triangle_2d(scalar x0, scalar y0, scalar x1, scalar y1,
                                 scalar x2, scalar y2, scalar x3, scalar y3,
                                 scalar &a, scalar &b, scalar &c) {
  x1 -= x0;
  x2 -= x0;
  x3 -= x0;
  y1 -= y0;
  y2 -= y0;
  y3 -= y0;
  int signa = orientation(x2, y2, x3, y3, a);
  if (signa == 0) return false;
  int signb = orientation(x3, y3, x1, y1, b);
  if (signb != signa) return false;
  int signc = orientation(x1, y1, x2, y2, c);
  if (signc != signa) return false;
  scalar sum = a + b + c;
  assert(sum != 0);  // if the SOS signs match and are nonkero, there's no way
                     // all of a, b, and c are zero.
  a /= sum;
  b /= sum;
  c /= sum;
  return true;
}

void make_level_set3(const std::vector<Vector3i> &tri,
                     const std::vector<Vector3s> &x, const Vector3s &origin,
                     scalar dx, int ni, int nj, int nk, Array3d &phi,
                     const int exact_band) {
  phi.resize(ni, nj, nk);
  phi.assign((ni + nj + nk) * dx);  // upper bound on distance
  Array3i closest_tri(ni, nj, nk, -1);
  Array3i intersection_count(ni, nj, nk,
                             0);  // intersection_count(i,j,k) is # of tri
                                  // intersections in (i-1,i]x{j}x{k}
  // we begin by initializing distances near the mesh, and figuring out
  // intersection counts
  Vector3s ijkmin, ijkmax;
  for (unsigned int t = 0; t < tri.size(); ++t) {
    const Vector3i &pqr = tri[t];
    // coordinates in grid to high precision
    scalar fip = (x[pqr(0)][0] - origin[0]) / dx,
           fjp = (x[pqr(0)][1] - origin[1]) / dx,
           fkp = (x[pqr(0)][2] - origin[2]) / dx;
    scalar fiq = (x[pqr(1)][0] - origin[0]) / dx,
           fjq = (x[pqr(1)][1] - origin[1]) / dx,
           fkq = (x[pqr(1)][2] - origin[2]) / dx;
    scalar fir = (x[pqr(2)][0] - origin[0]) / dx,
           fjr = (x[pqr(2)][1] - origin[1]) / dx,
           fkr = (x[pqr(2)][2] - origin[2]) / dx;
    // do distances nearby
    int i0 = mathutils::clamp(
        int(std::min(std::min(fip, fiq), fir)) - exact_band, 0, ni - 1);
    int i1 = mathutils::clamp(
        int(std::max(std::max(fip, fiq), fir)) + exact_band + 1, 0, ni - 1);
    int j0 = mathutils::clamp(
        int(std::min(std::min(fjp, fjq), fjr)) - exact_band, 0, nj - 1);
    int j1 = mathutils::clamp(
        int(std::max(std::max(fjp, fjq), fjr)) + exact_band + 1, 0, nj - 1);
    int k0 = mathutils::clamp(
        int(std::min(std::min(fkp, fkq), fkr)) - exact_band, 0, nk - 1);
    int k1 = mathutils::clamp(
        int(std::max(std::max(fkp, fkq), fkr)) + exact_band + 1, 0, nk - 1);
    for (int k = k0; k <= k1; ++k)
      for (int j = j0; j <= j1; ++j)
        for (int i = i0; i <= i1; ++i) {
          Vector3s gx(i * dx + origin[0], j * dx + origin[1],
                      k * dx + origin[2]);
          scalar d =
              point_triangle_distance(gx, x[pqr(0)], x[pqr(1)], x[pqr(2)]);
          if (d < phi(i, j, k)) {
            phi(i, j, k) = d;
            closest_tri(i, j, k) = t;
          }
        }
    // and do intersection counts
    j0 = mathutils::clamp((int)std::ceil(std::min(std::min(fjp, fjq), fjr)), 0,
                          nj - 1);
    j1 = mathutils::clamp((int)std::floor(std::max(std::max(fjp, fjq), fjr)), 0,
                          nj - 1);
    k0 = mathutils::clamp((int)std::ceil(std::min(std::min(fkp, fkq), fkr)), 0,
                          nk - 1);
    k1 = mathutils::clamp((int)std::floor(std::max(std::max(fkp, fkq), fkr)), 0,
                          nk - 1);
    for (int k = k0; k <= k1; ++k)
      for (int j = j0; j <= j1; ++j) {
        scalar a, b, c;
        if (point_in_triangle_2d(j, k, fjp, fkp, fjq, fkq, fjr, fkr, a, b, c)) {
          scalar fi = a * fip + b * fiq + c * fir;  // intersection i coordinate
          int i_interval = int(
              std::ceil(fi));  // intersection is in (i_interval-1,i_interval]
          if (i_interval < 0)
            ++intersection_count(0, j,
                                 k);  // we enlarge the first interval to
                                      // include everything to the -x direction
          else if (i_interval < ni)
            ++intersection_count(i_interval, j, k);
          // we ignore intersections that are beyond the +x side of the grid
        }
      }
  }
  // and now we fill in the rest of the distances with fast sweeping
  for (unsigned int pass = 0; pass < 2; ++pass) {
    sweep(tri, x, phi, closest_tri, origin, dx, +1, +1, +1);
    sweep(tri, x, phi, closest_tri, origin, dx, -1, -1, -1);
    sweep(tri, x, phi, closest_tri, origin, dx, +1, +1, -1);
    sweep(tri, x, phi, closest_tri, origin, dx, -1, -1, +1);
    sweep(tri, x, phi, closest_tri, origin, dx, +1, -1, +1);
    sweep(tri, x, phi, closest_tri, origin, dx, -1, +1, -1);
    sweep(tri, x, phi, closest_tri, origin, dx, +1, -1, -1);
    sweep(tri, x, phi, closest_tri, origin, dx, -1, +1, +1);
  }
  // then figure out signs (inside/outside) from intersection counts
  for (int k = 0; k < nk; ++k)
    for (int j = 0; j < nj; ++j) {
      int total_count = 0;
      for (int i = 0; i < ni; ++i) {
        total_count += intersection_count(i, j, k);
        if (total_count % 2 ==
            1) {  // if parity of intersections so far is odd,
          phi(i, j, k) = -phi(i, j, k);  // we are inside the mesh
        }
      }
    }
}
