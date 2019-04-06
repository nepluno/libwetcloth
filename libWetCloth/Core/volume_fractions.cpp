//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2008 Christopher Batty, Robert Bridson
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

#include "MathUtilities.h"
#include "volume_fractions.h"

// Assumes phi0<0 and phi1>=0, phi2>=0, or vice versa.
// In particular, phi0 must not equal either of phi1 or phi2.
template<class T>
static T sorted_triangle_fraction(T phi0, T phi1, T phi2)
{
  return phi0*phi0/(2*(phi0-phi1)*(phi0-phi2));
}

float area_fraction(float phi0, float phi1, float phi2)
{
  if(phi0<0){
    if(phi1<0){
      if(phi2<0)    return 0;
      else          return 1-sorted_triangle_fraction(phi2, phi0, phi1);
    }else if(phi2<0) return 1-sorted_triangle_fraction(phi1, phi2, phi0);
    else             return sorted_triangle_fraction(phi0, phi1, phi2);
  }else if(phi1<0){
    if(phi2<0)       return 1-sorted_triangle_fraction(phi0, phi1, phi2);
    else             return sorted_triangle_fraction(phi1, phi2, phi0);
  }else if(phi2<0)    return sorted_triangle_fraction(phi2, phi0, phi1);
  else                return 0;
}

double area_fraction(double phi0, double phi1, double phi2)
{
  if(phi0<0){
    if(phi1<0){
      if(phi2<0)    return 0;
      else          return 1-sorted_triangle_fraction(phi2, phi0, phi1);
    }else if(phi2<0) return 1-sorted_triangle_fraction(phi1, phi2, phi0);
    else             return sorted_triangle_fraction(phi0, phi1, phi2);
  }else if(phi1<0){
    if(phi2<0)       return 1-sorted_triangle_fraction(phi0, phi1, phi2);
    else             return sorted_triangle_fraction(phi1, phi2, phi0);
  }else if(phi2<0)    return sorted_triangle_fraction(phi2, phi0, phi1);
  else                return 0;
}

float area_fraction(float phi00, float phi10, float phi01, float phi11)
{
  float phimid=(phi00+phi10+phi01+phi11)/4;
  return (area_fraction(phi00, phi10, phimid)
          +area_fraction(phi10, phi11, phimid)
          +area_fraction(phi11, phi01, phimid)
          +area_fraction(phi01, phi00, phimid))/4;
}

double area_fraction(double phi00, double phi10, double phi01, double phi11)
{
  double phimid=(phi00+phi10+phi01+phi11)/4;
  return (area_fraction(phi00, phi10, phimid)
          +area_fraction(phi10, phi11, phimid)
          +area_fraction(phi11, phi01, phimid)
          +area_fraction(phi01, phi00, phimid))/4;
}

//============================================================================

// Assumes phi0<0 and phi1>=0, phi2>=0, and phi3>=0 or vice versa.
// In particular, phi0 must not equal any of phi1, phi2 or phi3.
template<class T>
static T sorted_tet_fraction(T phi0, T phi1, T phi2, T phi3)
{
  return phi0*phi0*phi0/((phi0-phi1)*(phi0-phi2)*(phi0-phi3));
}

// Assumes phi0<0, phi1<0, and phi2>=0, and phi3>=0 or vice versa.
// In particular, phi0 and phi1 must not equal any of phi2 and phi3.
template<class T>
static T sorted_prism_fraction(T phi0, T phi1, T phi2, T phi3)
{
  T a=phi0/(phi0-phi2),
  b=phi0/(phi0-phi3),
  c=phi1/(phi1-phi3),
  d=phi1/(phi1-phi2);
  return a*b*(1-d)+b*(1-c)*d+c*d;
}

float volume_fraction(float phi0, float phi1, float phi2, float phi3)
{
  mathutils::sort(phi0, phi1, phi2, phi3);
  if(phi3<=0) return 1;
  else if(phi2<=0) return 1-sorted_tet_fraction(phi3, phi2, phi1, phi0);
  else if(phi1<=0) return sorted_prism_fraction(phi0, phi1, phi2, phi3);
  else if(phi0<=0) return sorted_tet_fraction(phi0, phi1, phi2, phi3);
  else return 0;
}

double volume_fraction(double phi0, double phi1, double phi2, double phi3)
{
  mathutils::sort(phi0, phi1, phi2, phi3);
  if(phi3<=0) return 1;
  else if(phi2<=0) return 1-sorted_tet_fraction(phi3, phi2, phi1, phi0);
  else if(phi1<=0) return sorted_prism_fraction(phi0, phi1, phi2, phi3);
  else if(phi0<=0) return sorted_tet_fraction(phi0, phi1, phi2, phi3);
  else return 0;
}

float volume_fraction(float phi000, float phi100,
                      float phi010, float phi110,
                      float phi001, float phi101,
                      float phi011, float phi111)
{
  // This is the average of the two possible decompositions of the cube into
  // five tetrahedra.
  return (   volume_fraction(phi000, phi001, phi101, phi011)
          +volume_fraction(phi000, phi101, phi100, phi110)
          +volume_fraction(phi000, phi010, phi011, phi110)
          +volume_fraction(phi101, phi011, phi111, phi110)
          +2*volume_fraction(phi000, phi011, phi101, phi110)
          +volume_fraction(phi100, phi101, phi001, phi111)
          +volume_fraction(phi100, phi001, phi000, phi010)
          +volume_fraction(phi100, phi110, phi111, phi010)
          +volume_fraction(phi001, phi111, phi011, phi010)
          +2*volume_fraction(phi100, phi111, phi001, phi010))/12;
}

double volume_fraction(double phi000, double phi100,
                       double phi010, double phi110,
                       double phi001, double phi101,
                       double phi011, double phi111)
{
  // This is the average of the two possible decompositions of the cube into
  // five tetrahedra.
  return (   volume_fraction(phi000, phi001, phi101, phi011)
          +volume_fraction(phi000, phi101, phi100, phi110)
          +volume_fraction(phi000, phi010, phi011, phi110)
          +volume_fraction(phi101, phi011, phi111, phi110)
          +2*volume_fraction(phi000, phi011, phi101, phi110)
          +volume_fraction(phi100, phi101, phi001, phi111)
          +volume_fraction(phi100, phi001, phi000, phi010)
          +volume_fraction(phi100, phi110, phi111, phi010)
          +volume_fraction(phi001, phi111, phi011, phi010)
          +2*volume_fraction(phi100, phi111, phi001, phi010))/12;
}

