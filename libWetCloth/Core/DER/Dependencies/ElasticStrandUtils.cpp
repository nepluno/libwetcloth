//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "ElasticStrandUtils.h"
#include "../../TimingUtilities.h"
#include <numeric>
#include <cmath>

Vec3 parallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 )
{
    // Compute rotation axis (if any)
    Vec3 b = t0.cross( t1 );
    const scalar bNorm = b.norm();
    if( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3& n0 = t0.cross( b ).normalized();
    const Vec3& n1 = t1.cross( b ).normalized();

    return u.dot( t0.normalized() ) * t1.normalized() + u.dot( n0 ) * n1 + u.dot( b ) * b;
}

Vec3 normalParallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 )
{
    // This should be called only to transport an orthogonal vector
    assert( isSmall(u.dot(t0)) );

    // Compute rotation axis (if any)
    Vec3 b = t0.cross( t1 );
    const scalar bNorm = b.norm();
    if( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3& n0 = t0.cross( b ).normalized();
    const Vec3& n1 = t1.cross( b ).normalized();

    return u.dot( n0 ) * n1 + u.dot( b ) * b;
}

Vec3 orthonormalParallelTransport( const Vec3& u, const Vec3& t0, const Vec3& t1 )
{
    // This should be called only to transport an orthogonal vector
    assert( isSmall(u.dot(t0)) );

    Vec3 b = t0.cross( t1 );
    const scalar bNorm = b.norm();
    if( isSmall( bNorm ) ) // vectors are nearly collinear
        return u;
    b /= bNorm;

    const Vec3& n0 = t0.cross( b );
    const Vec3& n1 = t1.cross( b );

    return u.dot( n0 ) * n1 + u.dot( b ) * b;
}

bool containsNans( const VecX &dofs )
{
    bool panic = false;
    for( int i = 0; i < dofs.size(); ++i )
    {
        if( std::isnan( dofs[i] ) )
        {
            panic = true;
            break;
        }
    }

    return panic;
}


void genCurlyHair( const Vec3& initnorm, const Vec3& startpoint, const double& dL, const int& nv, std::vector<Vec3>& vertices, double curl_radius, double curl_density, double root_length )
{
  // generate an orthonormal frame
  Vec3 p1;
  p1 = findNormal(initnorm);
  Vec3 p2 = p1.cross( initnorm );
  
  vertices.push_back( startpoint );
  Vec3 freepoint( startpoint + root_length * initnorm );
  vertices.push_back( freepoint );
  // for each curve sample vertex
  scalar xa = M_PI/( curl_density * 4 ); // 0 // start curve parameter
  scalar xb = 0; // end curve parameter
  for( int j = 1; j < nv-1; ++j )
  {
    xb = ( dL + xa + curl_radius * curl_radius * curl_density * curl_density * xa)/(1 + curl_radius * curl_radius * curl_density * curl_density ); //upate to get length dL along curve
    vertices.push_back(
                       freepoint
                       + xb * initnorm
                       +  curl_radius * cos( xb * curl_density  ) * p2
                       +  curl_radius * sin( xb * curl_density  ) * p1
                       );
    xa = xb; // next...
  }
}

void updateCurlyHair( const double& dL, std::vector<Vec3>& vertices, double curl_radius, double curl_density, double root_length )
{
  const int nv = vertices.size();
  if(nv < 2) return;
  
  // generate an orthonormal frame
  Vec3 initnorm = (vertices[1] - vertices[0]).normalized();
  vertices[1] = vertices[0] + initnorm * root_length;
  
  Vec3 p1;
  p1 = findNormal(initnorm);
  Vec3 p2 = p1.cross( initnorm );
  
  scalar xa = M_PI/( curl_density * 4 ); // 0 // start curve parameter
  scalar xb = 0; // end curve parameter
  
  const Vec3& freepoint = vertices[1];
  
  for( int j = 2; j < nv; ++j )
  {
    xb = ( dL + xa + curl_radius * curl_radius * curl_density * curl_density * xa)/(1 + curl_radius * curl_radius * curl_density * curl_density ); //upate to get length dL along curve
    vertices[j] = freepoint
    + xb * initnorm
    +  curl_radius * cos( xb * curl_density  ) * p2
    +  curl_radius * sin( xb * curl_density  ) * p1;
    xa = xb; // next...
  }
}
