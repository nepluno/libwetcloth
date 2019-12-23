//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ELASTICSTRANDUTILS_H_
#define ELASTICSTRANDUTILS_H_

#include <limits>
#include <stdexcept>
#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../Definitions.h"

/**
 * \brief Tests if a matrix is symmetric
 */
template<typename MatrixT>
inline bool isSymmetric( const MatrixT& A )
{
    for( int i = 0; i < A.rows(); ++i ){
        for( int j = i + 1; j < A.cols(); ++j ){
            if( !isSmall( A( i, j ) - A( j, i ) ) )
            {
                std::cerr << "isSymmetric failing by " << fabs( A( i, j ) - A( j, i ) ) << '\n';
                return false;
            }
        }
    }
    return true;
}

/**
 * \brief Projects v on the plane normal to n and normalize it
 */
inline void orthoNormalize( Vec3& v, const Vec3& n )
{
    assert( isApproxUnit(n) );
    v -= v.dot( n ) * n;
    v.normalize();
}

/**
 * \brief Parallel-transports u along the t0->t1 transformation.
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t0 first vector defining transformation. It doesn't have to be normalized.
 * \param [in] t1 second vector defining transformation. It doesn't have to be normalized.
 *
 * \return If t0 and t1 are colinear, then u is unchanged. Otherwise (t0, t1) defines a plane
 * and a rotation around an axis perpendicular to that plane, which sends t0.normalized()
 * onto t1.normalized(). The parallel transport of u is its image by this rotation.
 *
 * \see normalParallelTransport()
 */
Vec3 parallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 );

/**
 * \brief Parallel-transports u along the t0->t1 transformation, assuming that u is normal to t0.
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t0 first vector defining transformation. It doesn't have to be normalized.
 * \param [in] t1 second vector defining transformation. It doesn't have to be normalized.
 *
 * \return If t0 and t1 are colinear, then u is unchanged. Otherwise (t0, t1) defines a plane
 * and a rotation around an axis perpendicular to that plane, which sends t0.normalized()
 * onto t1.normalized(). The parallel transport of u is its image by this rotation.
 *
 * \note It is assumed (see assertion below) that u.dot(0)==0, which is the case when
 * parallel-transporting frame vectors for instance. Then this is cheaper to call
 * than parallelTransport()
 */
Vec3 normalParallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 );

/**
 * \brief Parallel-transports u along the t0->t1 transformation, assuming that u is normal to t0.
 * and all the vectors are unitary
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t0 first vector defining transformation. It  have to be normalized.
 * \param [in] t1 second vector defining transformation. It have to be normalized.
 *
 * \return If t0 and t1 are colinear, then u is unchanged. Otherwise (t0, t1) defines a plane
 * and a rotation around an axis perpendicular to that plane, which sends t0
 * onto t1. The parallel transport of u is its image by this rotation.
 *
 * \note It is assumed that u.dot( t0 ) = 0 and t0.norm() = t1.norm() = u.norm() = 1
 */
Vec3 orthonormalParallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 );

void genCurlyHair( const Vec3& initnorm, const Vec3& startpoint, const double& dL, const int& nv, std::vector<Vec3>& vertices, double curl_radius, double curl_density, double root_length );

void updateCurlyHair( const double& dL, std::vector<Vec3>& vertices, double curl_radius, double curl_density, double root_length );
/**
 * \brief Finds an arbitrary unit vector orthogonal to u
 *
 * \param [in] u vector in dimension n
 * \return An arbitrary vector orthogonal to u
 */
inline Eigen::Matrix<scalar, Eigen::Dynamic, 1> findNormal( const Eigen::Matrix<scalar, Eigen::Dynamic, 1>& u )
{
    assert( u.norm() != 0 );
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> v(u.size());
	v.setZero();

    int maxCoordinate = 0;
	int n = u.size();
    for( int i = 0; i < n; ++i )
    {
        if( isSmall( u[i] ) )
        {
            v[i] = 1;
            goto finished;
        }
        if( fabs( u[i] ) > fabs( u[maxCoordinate] ) )
            maxCoordinate = i;
    }
    {
        const int otherCoordinate = ( maxCoordinate + 1 ) % n;
        v[otherCoordinate] = u[maxCoordinate];
        v[maxCoordinate] = -u[otherCoordinate];
    }
    v.normalize();

    finished: assert( isSmall( u.dot( v ) ) );

    return v;
}

/**
 * \brief Computes the signed angle from one vector to another.
 *
 * \param [in] u first input vector
 * \param [in] v second input vector
 * \param [in] n orientation vector
 *
 * The sign of the angle is positive if u.cross(v) is in the same half-space as n.
 */
inline scalar signedAngle( const Vec3& u, const Vec3& v, const Vec3& n )
{
    Vec3 w = u.cross( v );
    scalar angle = atan2( w.norm(), u.dot( v ) );
    if( n.dot( w ) < 0 )
        return -angle;
    return angle;
}

/**
 * \brief Rotates a vector
 *
 * \param [in] v vector to rotate
 * \param [in] z normalized vector on the rotation axis
 * \param [in] theta rotation angle
 *
 */
template<typename scalarT>
inline void rotateAxisAngle( typename Eigen::Matrix<scalarT, 3, 1> & v,
        const typename Eigen::Matrix<scalarT, 3, 1> & z, const scalarT theta )
{
    assert( isApproxUnit( z ) );

    if( theta == 0 )
        return;

    const scalarT c = cos( theta );
    const scalarT s = sin( theta );

    v = c * v + s * z.cross( v ) + z.dot( v ) * ( 1.0 - c ) * z;
}

/**
 * \brief Outer product of two vectors.
 */
template<int n>
inline Eigen::Matrix<scalar, n, n> outerProd( const Eigen::Matrix<scalar, n, 1>& a,
        const Eigen::Matrix<scalar, n, 1>& b )
{
    return a * b.transpose();
}

/**
 * \brief Matrix representation of the cross product operator.
 */
inline Mat3 crossMat( const Vec3& a )
{
    Mat3 M;
    M << 0, -a( 2 ), a( 1 ), a( 2 ), 0, -a( 0 ), -a( 1 ), a( 0 ), 0;

    return M;
}

/**
 * \brief Computes u^T B v, assuming B is symmetric 2x2 and u, v are 2x1 vectors.
 */
inline scalar innerBProduct( const Mat2& B, const Vec2& u, const Vec2& v )
{
    assert( isSymmetric( B ) );

    return B( 0, 0 ) * u[0] * v[0] + B( 0, 1 ) * ( u[0] * v[1] + u[1] * v[0] )
            + B( 1, 1 ) * u[1] * v[1]; // Good
}

/**
 * \brief Computes Q B Q^T, assuming B is symmetric 2x2 and Q is nx2. The result is then (exactly) symmetric nxn.
 */
template<int n>
inline void symBProduct( Eigen::Matrix<scalar, n, n>& result, const Mat2& B,
        const Eigen::Matrix<scalar, n, 2>& Q )
{
    assert( isSymmetric( B ) );

    for( int i = 0; i < n; ++i )
    {
        const Vec2& Qrow_i = Q.row( i );
        result( i, i ) = innerBProduct( B, Qrow_i, Qrow_i );
        for( int j = 0; j < i; ++j )
            result( i, j ) = result( j, i ) = innerBProduct( B, Qrow_i, Q.row( j ) );
    }
}

/**
 * Angular interpolation by a factor t between two vectors v0 and v1
 */
template<typename VectorT, typename scalarT>
inline VectorT vectorSlerp( const VectorT& v0, const VectorT &v1, scalarT t )
{
    const scalarT angle = std::acos( clamp( v0.dot( v1 ), -1., 1. ) );
    if( isSmall( angle ) )
        return v0;
    const scalarT invSin = 1. / std::sin( angle );
    return invSin * ( std::sin( ( 1. - t ) * angle ) * v0 + std::sin( t * angle ) * v1 );
}

/**
 * Applies the (projective space) transform M to x identified to (x,1)
 */
template<typename scalarT>
Eigen::Matrix<scalarT, 3, 1> transformPoint( const Eigen::Matrix<scalarT, 4, 4>& M,
        const Eigen::Matrix<scalarT, 3, 1> & x )
{
    return M.template block<3, 3>( 0, 0 ) * x + M.template block<3, 1>( 0, 3 );
}

/**
 * Tests whether dofs contains NaNs
 */
bool containsNans( const VecX &dofs );

/**
 \return an angle theta in [-pi, pi] such that theta == \p angle mod (2pi)
 */
inline scalar clamp2Pi( scalar angle )
{
    scalar theta = angle ;
    while ( theta > M_PI )
        theta -= 2. * M_PI;
    while ( theta <= -M_PI )
        theta += 2. * M_PI;

    return theta ;
}

#endif
