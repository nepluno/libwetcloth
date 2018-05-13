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

#ifndef CAMERA_H
#define CAMERA_H
#include <Eigen/Dense>
#include "MathDefs.h"

class Camera
{
public:
  Eigen::Quaterniond rotation_;	// rotation
  Eigen::Vector3d center_;	// center point
  double	 dist_;			// point of view to center
  double   radius_;		// bounding sphere
  double   fov_;			// angle
  Camera ( const Camera& that );
  void operator = ( const Camera& that );

  explicit Camera ( const double fov = 40 );
  explicit Camera ( Eigen::Quaterniond& rot, Eigen::Vector3d& center, const double& dist, const double& radius, const double& fov );
  void init( const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax );
  void clone ( const Camera& that );
  void getViewDir( Eigen::Vector3d& viewdir ) const;
  void getLookAt( Eigen::Vector3d& eye, Eigen::Vector3d& center, Eigen::Vector3d& up ) const;
  void getEye( Eigen::Vector3d& eye ) const;
  void getPerspective( double& fov, double& zNear, double& zFar ) const;
  void rotate( const double  oldx, const double oldy, const double newx, const double newy );
  void zoom( const double  oldx, const double oldy, const double newx, const double newy );
  void pan( const double  oldx, const double oldy, const double newx, const double newy );
  void rotate( const Vector3s& axis, const scalar& angle, bool global );
  void project_to_sphere( const double& radius, Eigen::Vector3d& p ) const;
  friend std::ostream &operator<<( std::ostream &output, const Camera &cam );
};
#endif//CAMERA_H
