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

#ifndef __RENDERING_UTILITIES_H__
#define __RENDERING_UTILITIES_H__

#ifdef WIN32
#include <Windows.h>
#endif

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <list>
#include <iostream>
#include <cstdio>

#include "MathDefs.h"
#include "StringUtilities.h"
#include "MathUtilities.h"
#include "JET.h"

namespace renderingutils
{
	// False => error
	bool checkGLErrors();
	
	class Color
	{
	public:
		
		Color();
		
		Color( double r, double g, double b );
		
		Color( const Vector3s& );
		
		Vector3s toVector() const;
		
		double r;
		double g;
		double b;
	};
	
	struct Viewport
	{
	public:
		Viewport() : cx(0.0), cy(0.0), cz(0.0), size(0.) {}
		double cx;
		double cy;
		double cz;
		double rx;
		double ry;
		double rz;
		double size;
	};
	
	inline Vector3s interpolateColor(const scalar& x, const scalar xmin = 0.0, const scalar xmax = 1.0)
	{
		scalar dm = (xmax - xmin);
		
		scalar a;
		if(dm == 0.0) a = x;
		else a = (x - xmin) / dm * (scalar)(jetmapping_size - 1);
		
		int isel = std::max(std::min((int) a, jetmapping_size - 1), 0);
		int inext = (isel + 1) % (jetmapping_size);
		scalar fraca = std::max(std::min(a - (scalar) isel, 1.0), 0.0);
		
		return mathutils::lerp(jetmapping_real[isel], jetmapping_real[inext], fraca);
	}
}

#endif
