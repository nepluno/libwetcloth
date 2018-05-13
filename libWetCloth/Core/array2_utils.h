//
// This file is part of the libWetCloth open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// We would like to hear from you if you appreciate this work.
//
// The Clear BSD License
//
// Copyright 2008 Christopher Batty, Robert Bridson
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
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


#ifndef ARRAY2_UTILS_H
#define ARRAY2_UTILS_H

#include "array2.h"
#include "MathUtilities.h"

using namespace mathutils;

template<class S, class T>
T interpolate_value(const Eigen::Matrix<S, 2, 1>& point, const Array2<T, Array1<T> >& grid) {
   int i,j;
   S fx,fy;

   get_barycentric(point[0], i, fx, 0, grid.ni);
   get_barycentric(point[1], j, fy, 0, grid.nj);
   
   return bilerp(
      grid(i,j), grid(i+1,j), 
      grid(i,j+1), grid(i+1,j+1), 
      fx, fy);
}

template<class S, class T>
T interpolate_value(const Eigen::Matrix<S, 2, 1>& point, const std::vector<T>& grid, int ni, int nj) {
  int i,j;
  S fx,fy;
  
  get_barycentric(point[0], i, fx, 0, ni);
  get_barycentric(point[1], j, fy, 0, nj);
  
  return bilerp(
                grid[j * ni + i], grid[j * ni + i + 1],
                grid[(j + 1) * ni + i], grid[(j + 1) * ni + (i + 1)],
                fx, fy);
}


template<class T>
Eigen::Matrix<T, 2, 1> affine_interpolate_value(const Eigen::Matrix<T, 2, 1>& point, const Array2<T, Array1<T> >& grid) {
  int i,j;
  T fx,fy;
  
  get_barycentric(point[0], i, fx, 0, grid.ni);
  get_barycentric(point[1], j, fy, 0, grid.nj);
  
  return grad_bilerp(
                     grid(i,j), grid(i+1,j),
                     grid(i,j+1), grid(i+1,j+1),
                     fx, fy);
}

template<class S, class T>
T interpolate_gradient(Eigen::Matrix<T, 2, 1>& gradient, const Eigen::Matrix<S, 2, 1>& point, const Array2<T, Array1<T> >& grid) {
   int i,j;
   S fx,fy;
   get_barycentric(point[0], i, fx, 0, grid.ni);
   get_barycentric(point[1], j, fy, 0, grid.nj);

   T v00 = grid(i,j);
   T v01 = grid(i,j+1);
   T v10 = grid(i+1,j);
   T v11 = grid(i+1,j+1);

   T ddy0 = (v01 - v00);
   T ddy1 = (v11 - v10);

   T ddx0 = (v10 - v00);
   T ddx1 = (v11 - v01);

   gradient[0] = lerp(ddx0,ddx1,fy);
   gradient[1] = lerp(ddy0,ddy1,fx);
   
   //may as well return value too
   return bilerp(v00, v10, v01, v11, fx, fy);
}

template<typename S, typename T>
void interpolate_gradient(Eigen::Matrix<S, 2, 1>& gradient, const Eigen::Matrix<S, 2, 1>& point, const Array2<T, Array1<T> >& grid, S cellSize)
{
  int i,j;
  S fx,fy;
  get_barycentric(point(0), i, fx, 0, grid.ni);
  get_barycentric(point(1), j, fy, 0, grid.nj);
  
  S v00 = (S) grid(i,j);
  S v01 = (S) grid(i,j+1);
  S v10 = (S) grid(i+1,j);
  S v11 = (S) grid(i+1,j+1);
  
  S ddy0 = (v01 - v00);
  S ddy1 = (v11 - v10);
  
  S ddx0 = (v10 - v00);
  S ddx1 = (v11 - v01);
  
  gradient[0] = lerp(ddx0,ddx1,fy) / cellSize;
  gradient[1] = lerp(ddy0,ddy1,fx) / cellSize;
}

template<typename S, typename T>
void interpolate_gradient(Eigen::Matrix<S, 2, 1>& gradient, const Eigen::Matrix<S, 2, 1>& point, const std::vector<T>& grid, int ni, int nj, S cellSize)
{
  int i,j;
  S fx,fy;
  get_barycentric(point(0), i, fx, 0, ni);
  get_barycentric(point(1), j, fy, 0, nj);
  
  S v00 = (S) grid[j * ni + i];
  S v01 = (S) grid[(j + 1) * ni + i];
  S v10 = (S) grid[j * ni + i + 1];
  S v11 = (S) grid[(j + 1) * ni + i + 1];
  
  S ddy0 = (v01 - v00);
  S ddy1 = (v11 - v10);
  
  S ddx0 = (v10 - v00);
  S ddx1 = (v11 - v01);
  
  gradient[0] = lerp(ddx0,ddx1,fy) / cellSize;
  gradient[1] = lerp(ddy0,ddy1,fx) / cellSize;
}

template<typename S, typename T>
void interpolate_hessian(Eigen::Matrix<S, 2, 2>& hessian, const Eigen::Matrix<S, 2, 1>& point, const Array2<T, Array1<T> >& grid, S cellSize) {
  Eigen::Matrix<S, 2, 1> gx0;
  Eigen::Matrix<S, 2, 1> gx1;
  Eigen::Matrix<S, 2, 1> gy0;
  Eigen::Matrix<S, 2, 1> gy1;
  
  Eigen::Matrix<S, 2, 1> px0 = point - Eigen::Matrix<S, 2, 1>(cellSize * (S) 0.5, (S) 0.0);
  Eigen::Matrix<S, 2, 1> px1 = point + Eigen::Matrix<S, 2, 1>(cellSize * (S) 0.5, (S) 0.0);
  Eigen::Matrix<S, 2, 1> py0 = point - Eigen::Matrix<S, 2, 1>((S) 0.0, cellSize * (S) 0.5);
  Eigen::Matrix<S, 2, 1> py1 = point + Eigen::Matrix<S, 2, 1>((S) 0.0, cellSize * (S) 0.5);
  
  
  interpolate_gradient(gx0, px0, grid, cellSize);
  interpolate_gradient(gx1, px1, grid, cellSize);
  interpolate_gradient(gy0, py0, grid, cellSize);
  interpolate_gradient(gy1, py1, grid, cellSize);
  
  hessian.col(0) = (gx1 - gx0) / cellSize;
  hessian.col(1) = (gy1 - gy0) / cellSize;
  hessian = (hessian + hessian.transpose()) * 0.5;
}

template<class T>
void write_matlab_array(std::ostream &output, Array2<T, Array1<T> >&a, const char *variable_name, bool transpose=false)
{
   output<<variable_name<<"=[";
   for(int j = 0; j < a.nj; ++j) {
      for(int i = 0; i < a.ni; ++i)  {
         output<<a(i,j)<<" ";
      }
      output<<";";
   }
   output<<"]";
   if(transpose)
      output<<"'";
   output<<";"<<std::endl;
}

#endif
