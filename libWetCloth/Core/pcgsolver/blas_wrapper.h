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

#ifndef BLAS_WRAPPER_H
#define BLAS_WRAPPER_H

// Simple placeholder code for BLAS calls - replace with calls to a real BLAS library

#include <vector>
#include <Eigen/Core>

namespace robertbridson {

namespace BLAS{
// dot products ==============================================================

inline double dot(const std::vector<double> &x, const std::vector<double> &y)
{
   double sum = 0;
   size_t n = x.size() < y.size() ? x.size() : y.size();
   return Eigen::Map<Eigen::VectorXd>((double*) &y[0], n).dot(Eigen::Map<Eigen::VectorXd>((double*) &x[0], n));
}

// inf-norm (maximum absolute value: index of max returned) ==================

inline int index_abs_max(const std::vector<double> &x)
{
   int maxind = 0;
   int col = 0;

   size_t n = x.size();
   Eigen::Map<Eigen::VectorXd>((double*) &x[0], n).cwiseAbs().maxCoeff(&maxind, &col);

   return maxind;
}

// inf-norm (maximum absolute value) =========================================
// technically not part of BLAS, but useful

inline double abs_max(const std::vector<double> &x)
{ return std::fabs(x[index_abs_max(x)]); }

// saxpy (y=alpha*x+y) =======================================================

inline void add_scaled(double alpha, const std::vector<double> &x, std::vector<double> &y)
{
    size_t n = x.size() < y.size() ? x.size() : y.size();
    Eigen::Map<Eigen::VectorXd>((double*) &y[0], n) += Eigen::Map<const Eigen::VectorXd>((const double*) &x[0], n) * alpha;
}
}
}
#endif
