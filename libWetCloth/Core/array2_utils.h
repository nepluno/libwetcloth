//
// This file is part of the libWetCloth open source project
//
// Copyright 2008 Christopher Batty, Robert Bridson
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ARRAY2_UTILS_H
#define ARRAY2_UTILS_H

#include "MathUtilities.h"
#include "array2.h"

using namespace mathutils;

template <class S, class T>
T interpolate_value(const Eigen::Matrix<S, 2, 1>& point,
                    const Array2<T, Array1<T> >& grid) {
  int i, j;
  S fx, fy;

  get_barycentric(point[0], i, fx, 0, grid.ni);
  get_barycentric(point[1], j, fy, 0, grid.nj);

  return bilerp(grid(i, j), grid(i + 1, j), grid(i, j + 1), grid(i + 1, j + 1),
                fx, fy);
}

template <class S, class T>
T interpolate_value(const Eigen::Matrix<S, 2, 1>& point,
                    const std::vector<T>& grid, int ni, int nj) {
  int i, j;
  S fx, fy;

  get_barycentric(point[0], i, fx, 0, ni);
  get_barycentric(point[1], j, fy, 0, nj);

  return bilerp(grid[j * ni + i], grid[j * ni + i + 1], grid[(j + 1) * ni + i],
                grid[(j + 1) * ni + (i + 1)], fx, fy);
}

template <class T>
Eigen::Matrix<T, 2, 1> affine_interpolate_value(
    const Eigen::Matrix<T, 2, 1>& point, const Array2<T, Array1<T> >& grid) {
  int i, j;
  T fx, fy;

  get_barycentric(point[0], i, fx, 0, grid.ni);
  get_barycentric(point[1], j, fy, 0, grid.nj);

  return grad_bilerp(grid(i, j), grid(i + 1, j), grid(i, j + 1),
                     grid(i + 1, j + 1), fx, fy);
}

template <class S, class T>
T interpolate_gradient(Eigen::Matrix<T, 2, 1>& gradient,
                       const Eigen::Matrix<S, 2, 1>& point,
                       const Array2<T, Array1<T> >& grid) {
  int i, j;
  S fx, fy;
  get_barycentric(point[0], i, fx, 0, grid.ni);
  get_barycentric(point[1], j, fy, 0, grid.nj);

  T v00 = grid(i, j);
  T v01 = grid(i, j + 1);
  T v10 = grid(i + 1, j);
  T v11 = grid(i + 1, j + 1);

  T ddy0 = (v01 - v00);
  T ddy1 = (v11 - v10);

  T ddx0 = (v10 - v00);
  T ddx1 = (v11 - v01);

  gradient[0] = lerp(ddx0, ddx1, fy);
  gradient[1] = lerp(ddy0, ddy1, fx);

  // may as well return value too
  return bilerp(v00, v10, v01, v11, fx, fy);
}

template <typename S, typename T>
void interpolate_gradient(Eigen::Matrix<S, 2, 1>& gradient,
                          const Eigen::Matrix<S, 2, 1>& point,
                          const Array2<T, Array1<T> >& grid, S cellSize) {
  int i, j;
  S fx, fy;
  get_barycentric(point(0), i, fx, 0, grid.ni);
  get_barycentric(point(1), j, fy, 0, grid.nj);

  S v00 = (S)grid(i, j);
  S v01 = (S)grid(i, j + 1);
  S v10 = (S)grid(i + 1, j);
  S v11 = (S)grid(i + 1, j + 1);

  S ddy0 = (v01 - v00);
  S ddy1 = (v11 - v10);

  S ddx0 = (v10 - v00);
  S ddx1 = (v11 - v01);

  gradient[0] = lerp(ddx0, ddx1, fy) / cellSize;
  gradient[1] = lerp(ddy0, ddy1, fx) / cellSize;
}

template <typename S, typename T>
void interpolate_gradient(Eigen::Matrix<S, 2, 1>& gradient,
                          const Eigen::Matrix<S, 2, 1>& point,
                          const std::vector<T>& grid, int ni, int nj,
                          S cellSize) {
  int i, j;
  S fx, fy;
  get_barycentric(point(0), i, fx, 0, ni);
  get_barycentric(point(1), j, fy, 0, nj);

  S v00 = (S)grid[j * ni + i];
  S v01 = (S)grid[(j + 1) * ni + i];
  S v10 = (S)grid[j * ni + i + 1];
  S v11 = (S)grid[(j + 1) * ni + i + 1];

  S ddy0 = (v01 - v00);
  S ddy1 = (v11 - v10);

  S ddx0 = (v10 - v00);
  S ddx1 = (v11 - v01);

  gradient[0] = lerp(ddx0, ddx1, fy) / cellSize;
  gradient[1] = lerp(ddy0, ddy1, fx) / cellSize;
}

template <typename S, typename T>
void interpolate_hessian(Eigen::Matrix<S, 2, 2>& hessian,
                         const Eigen::Matrix<S, 2, 1>& point,
                         const Array2<T, Array1<T> >& grid, S cellSize) {
  Eigen::Matrix<S, 2, 1> gx0;
  Eigen::Matrix<S, 2, 1> gx1;
  Eigen::Matrix<S, 2, 1> gy0;
  Eigen::Matrix<S, 2, 1> gy1;

  Eigen::Matrix<S, 2, 1> px0 =
      point - Eigen::Matrix<S, 2, 1>(cellSize * (S)0.5, (S)0.0);
  Eigen::Matrix<S, 2, 1> px1 =
      point + Eigen::Matrix<S, 2, 1>(cellSize * (S)0.5, (S)0.0);
  Eigen::Matrix<S, 2, 1> py0 =
      point - Eigen::Matrix<S, 2, 1>((S)0.0, cellSize * (S)0.5);
  Eigen::Matrix<S, 2, 1> py1 =
      point + Eigen::Matrix<S, 2, 1>((S)0.0, cellSize * (S)0.5);

  interpolate_gradient(gx0, px0, grid, cellSize);
  interpolate_gradient(gx1, px1, grid, cellSize);
  interpolate_gradient(gy0, py0, grid, cellSize);
  interpolate_gradient(gy1, py1, grid, cellSize);

  hessian.col(0) = (gx1 - gx0) / cellSize;
  hessian.col(1) = (gy1 - gy0) / cellSize;
  hessian = (hessian + hessian.transpose()) * 0.5;
}

template <class T>
void write_matlab_array(std::ostream& output, Array2<T, Array1<T> >& a,
                        const char* variable_name, bool transpose = false) {
  output << variable_name << "=[";
  for (int j = 0; j < a.nj; ++j) {
    for (int i = 0; i < a.ni; ++i) {
      output << a(i, j) << " ";
    }
    output << ";";
  }
  output << "]";
  if (transpose) output << "'";
  output << ";" << std::endl;
}

#endif
