//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2012 Jean-Marie Aubry
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
