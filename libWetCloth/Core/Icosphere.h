//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ICOSPHERE_H
#define ICOSPHERE_H

/* Icosphere generation algorithm.
 * Adapted from Andreas Kahler's C# implementation found here:
 * http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
 */

#include <unordered_map>
#include <vector>

#include "MathDefs.h"
#include "SolidMesh.h"

class Icosphere : public SolidMesh {
  std::unordered_map<uint64, int> middlePointIndexCache;
  int index;

  int addVertex(const Vector3s& p);
  int getMiddlePoint(int p1, int p2);

 public:
  Icosphere(int recursionLevel, const scalar& radius);
};

#endif
