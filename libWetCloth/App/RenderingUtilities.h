//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef RENDERING_UTILITIES_H
#define RENDERING_UTILITIES_H

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
	if (dm == 0.0) a = x;
	else a = (x - xmin) / dm * (scalar)(jetmapping_size - 1);

	int isel = std::max(std::min((int) a, jetmapping_size - 1), 0);
	int inext = (isel + 1) % (jetmapping_size);
	scalar fraca = std::max(std::min(a - (scalar) isel, 1.0), 0.0);

	return mathutils::lerp(jetmapping_real[isel], jetmapping_real[inext], fraca);
}
}

#endif
