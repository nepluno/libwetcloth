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
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
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

#ifndef _DER_DEFINITIONS_H_
#define _DER_DEFINITIONS_H_

#include <iostream>

#ifndef INCLUDE_VANILLA_EIGEN
#define EIGEN_VECTOR_IO_FORMAT Eigen::IOFormat(8, Eigen::DontAlignCols, ", ", ", ", "", "", "{ ", " }")
#define EIGEN_MATRIX_IO_FORMAT Eigen::IOFormat(8, 0, ", ", "\n", "{ ", " }", "{ ", " }")
#undef EIGEN_DEFAULT_IO_FORMAT // < To silence some warnings about redefining
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_VECTOR_IO_FORMAT
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#undef EIGEN_INITIALIZE_MATRICES_BY_ZERO // < To silence some warnings about redefining
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#endif // INCLUDE_VANILLA_EIGEN

#undef Success // Conflicts with Eigen
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

namespace Eigen
{
    template<typename _scalar, int _Options, typename _Index>
    class SparseMatrix;
}

typedef double scalar; ///< the scalar type
typedef uint16_t IndexType; ///< large unsigned int for IDs

typedef Eigen::Matrix<scalar, 2, 1> Vec2; ///< 2d scalar vector
typedef Eigen::Matrix<scalar, 3, 1> Vec3; ///< 3d scalar vector
typedef Eigen::Matrix<scalar, 4, 1> Vec4; ///< 3d scalar vector
typedef Eigen::Matrix<scalar, 11, 1> Vec11; ///< 11d scalar vector (stencil for local forces)
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VecX;

typedef std::vector<Vec2, Eigen::aligned_allocator<Vec2> > Vec2Array; ///< an array of 2d scalar vectors
typedef std::vector<Vec3, Eigen::aligned_allocator<Vec3> > Vec3Array;
typedef std::vector<Vec11, Eigen::aligned_allocator<Vec11> > Vec11Array; ///< an array of 11d scalar vectors

typedef Eigen::Matrix<scalar, 2, 2> Mat2; ///< 2x2 scalar matrix
typedef Eigen::Matrix<scalar, 3, 3> Mat3; ///< 3x3 scalar matrix
typedef Eigen::Matrix<scalar, 11, 11> Mat11; ///< 11x11 scalar matrix (stencil for local forces)
typedef std::pair<Mat11, Mat11> Mat11Pair;
typedef std::vector<Mat11, Eigen::aligned_allocator<Mat11> > Mat11Array; ///< an array of 11d scalar matrices

typedef Eigen::Triplet<scalar> Triplets;
typedef std::vector< Triplets > TripletXs;

template<typename scalarT>
scalarT EIGEN_STRONG_INLINE SMALL_NUMBER()
{
    return std::numeric_limits< scalarT >::epsilon() ;
}

template<>
float EIGEN_STRONG_INLINE SMALL_NUMBER<float>()
{
    return 1e-6;
}

template<>
double EIGEN_STRONG_INLINE SMALL_NUMBER<double>()
{
    return 1e-12;
}

EIGEN_STRONG_INLINE scalar square( const scalar x )
{
    return x * x;
}

EIGEN_STRONG_INLINE scalar cube( const scalar x )
{
    return x * x * x;
}

template<typename ComparableT>
EIGEN_STRONG_INLINE ComparableT clamp( const ComparableT x, const ComparableT l, const ComparableT u )
{
    return ( x > u ) ? u : ( ( x > l ) ? x : l );
}

template<typename scalarT>
EIGEN_STRONG_INLINE bool isSmall( scalarT x )
{
  return fabs( x ) < SMALL_NUMBER<scalarT>();
}

template<typename NormableT>
EIGEN_STRONG_INLINE bool isClose( const NormableT& x1, const NormableT& x2 )
{
    return isSmall( ( x1 - x2 ).norm() );
}

template<typename NormableT>
EIGEN_STRONG_INLINE bool isApproxUnit( const NormableT& x )
{
    return isSmall( x.squaredNorm() - 1 );
}

namespace std
{
    template < typename Derived >
    inline void swap ( Eigen::DenseBase< Derived >& a, Eigen::DenseBase< Derived >& b )
    {
        a.swap( b ) ;
    }

    template < typename Derived >
    inline void swap ( pair< Eigen::DenseBase< Derived >, Eigen::DenseBase< Derived > >& a,
                       pair< Eigen::DenseBase< Derived >, Eigen::DenseBase< Derived > >& b )
    {
        a.first.swap( b.first ) ;
        a.second.swap( b.second ) ;
    }
}

#endif 
