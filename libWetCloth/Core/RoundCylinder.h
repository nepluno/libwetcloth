//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ROUND_CYLINDER_H
#define ROUND_CYLINDER_H

/* Rounded Cylinder generation algorithm.
 */

#include "SolidMesh.h"
#include <vector>
#include <unordered_map>

class RoundCylinder : public SolidMesh
{
public:
	RoundCylinder(int N, int M, const double& ra, const double& rb, const double& h);
};



#endif
