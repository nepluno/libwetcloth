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


#ifndef VOLUME_FRACTIONS_H
#define VOLUME_FRACTIONS_H

// Given a triangle with level set values, use linear interpolation to
// estimate the fraction of the triangle occupied by the phi<0 part
float area_fraction(float phi0, float phi1, float phi2);
double area_fraction(double phi0, double phi1, double phi2);

// Given a rectangle with level set values, estimate fraction occupied
// by the phi<0 part
float area_fraction(float phi00, float phi10, float phi01, float phi11);
double area_fraction(double phi00, double phi10, double phi01, double phi11);

// Given a tetrahedron with level set values, use linear interpolation to
// estimate the fraction of the tetrahedron occupied by the phi<0 part
float volume_fraction(float phi0, float phi1, float phi2, float phi3);
double volume_fraction(double phi0, double phi1, double phi2, double phi3);

// Given a parallelepiped (e.g. cube) with level set values, estimate
// fraction occupied by the phi<0 part
float volume_fraction(float phi000, float phi100,
                      float phi010, float phi110,
                      float phi001, float phi101,
                      float phi011, float phi111);
double volume_fraction(double phi000, double phi100,
                       double phi010, double phi110,
                       double phi001, double phi101,
                       double phi011, double phi111);

#endif
