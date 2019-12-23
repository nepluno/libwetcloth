//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Camera.h"
#include <iostream>

Camera::Camera ( const double fov ): rotation_ ( Eigen::Quaterniond( 1,0,0,0 ) ), center_( Eigen::Vector3d() ), dist_( 0 ), radius_( 100 ), fov_( fov )
{
  Eigen::Vector3d b( 1,1,1 );
  b.normalize();
  b*= 100;
  this->init( -b, b );
}

Camera::Camera( const Camera& that )
{
  clone(that);
}

Camera::Camera ( Eigen::Quaterniond& rot, Eigen::Vector3d& center, const double& dist, const double& radius, const double& fov )
: rotation_(rot), center_(center), dist_(dist), radius_(radius), fov_(fov)
{}

void
Camera::init( const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax )
{
  this->center_ = ( bmin + bmax ) * 0.5;
  this->radius_ = ( bmin - bmax ).norm()* 0.5;
  this->dist_ = this->radius_  / std::sin ( this->fov_ / 360.0 * 3.1415926 );
  
  return;
}
void
Camera::clone ( const Camera& that )
{
  this->rotation_ = that.rotation_;
  this->center_ = that.center_;
  this->dist_   = that.dist_;
  this->radius_ = that.radius_;
  this->fov_    = that.fov_;
  return;
}

void Camera::operator = ( const Camera& that )
{
  clone(that);
}

void Camera::getViewDir( Eigen::Vector3d& viewdir ) const
{
  Eigen::Matrix3d r; // rotation matrix
  r.setIdentity();
  r = this->rotation_.toRotationMatrix();
  Eigen::Vector3d dir( 0, 0, 1.0 );
  viewdir = r * dir;
}

void
Camera::getLookAt( Eigen::Vector3d& eye, Eigen::Vector3d& center, Eigen::Vector3d& up ) const
{
  Eigen::Matrix3d r; // rotation matrix
  r.setIdentity();
  r = this->rotation_.toRotationMatrix();
  Eigen::Vector3d eye0( 0, 0, this->dist_ );
  Eigen::Vector3d up0( 0, 1, 0 );
  center = this->center_;
  
  eye = r * eye0 + this->center_;
  up  = r * up0;
  return;
}

void
Camera::getEye( Eigen::Vector3d& eye ) const
{
  Eigen::Matrix3d r; // rotation matrix
  r.setIdentity();
  r = this->rotation_.toRotationMatrix();
  Eigen::Vector3d eye0( 0, 0, this->dist_ );
  eye = r * eye0 + this->center_;
}

void
Camera::getPerspective( double& fov, double& zNear, double& zFar ) const
{
  fov = this->fov_;
  zNear = this->dist_ - this->radius_ * 10.0;
  zFar  = this->dist_ + this->radius_ * 10.0;
  if( zNear < 0 ) zNear = 0.1;
  
  return;
}
void
Camera::rotate( const double  oldx, const double oldy, const double newx, const double newy )
{
  Eigen::Vector3d oldp( oldx, oldy, 0 );
  Eigen::Vector3d newp( newx, newy, 0 );
  
  if ( oldp.isApprox( newp, 1.0e-16 ) ) return;
  
  double radius_virtual_sphere = 0.9;
  this->project_to_sphere( radius_virtual_sphere, oldp );
  this->project_to_sphere( radius_virtual_sphere, newp );
  Eigen::Quaterniond dr;
  dr.setFromTwoVectors( newp, oldp );
  this->rotation_ *= dr;
  return;
}

void
Camera::rotate( const Vector3s& axis, const scalar& angle, bool global )
{
  if(global)
  {
    this->rotation_ = Eigen::AngleAxis<scalar>(angle, axis) * this->rotation_;
  } else {
    this->rotation_ = this->rotation_ * Eigen::AngleAxis<scalar>(angle, axis);
  }
}

void
Camera::project_to_sphere( const double& radius, Eigen::Vector3d& p ) const
{
  p.z() = 0;
  const double d = p.x()* p.x()+ p.y() * p.y();
  const double r = radius * radius;
  if ( d < r )	p.z() = std::sqrt( r - d );
  else		p *= radius / p.norm();
  return;
}

void Camera::zoom( const double  oldx, const double oldy, const double newx, const double newy )
{
  this->radius_ = this->radius_ * (1.0 - newy + oldy);
  this->dist_ = this->radius_  / std::sin ( this->fov_ / 360.0 * 3.1415926 );
}

void Camera::pan( const double  oldx, const double oldy, const double newx, const double newy )
{
  Eigen::Matrix3d r; // rotation matrix
  r.setIdentity();
  r = this->rotation_.toRotationMatrix();
  
  Eigen::Vector3d up0( 0, 1, 0 );
  Eigen::Vector3d right0( 1, 0, 0);
  
  center_ -= r * (up0 * (newy - oldy) + right0 * (newx - oldx));
}

std::ostream& operator<<( std::ostream &output, const Camera &cam )
{
  output << "  <camera dist=\"" << cam.dist_ << "\" radius=\"" << cam.radius_ << "\" fov=\"" << cam.fov_ << "\">\n";
  output << "    <rotation x=\"" << cam.rotation_.x() << "\" y=\"" << cam.rotation_.y() << "\" z=\"" << cam.rotation_.z() << "\" w=\"" << cam.rotation_.w() << "\"/>\n";
  output << "    <center x=\"" << cam.center_.x() << "\" y=\"" << cam.center_.y() << "\" z=\"" << cam.center_.z() << "\"/>\n";
  output << "  </camera>";
  
  return output;
}
