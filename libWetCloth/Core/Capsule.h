//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef CAPSULE_H
#define CAPSULE_H

/* Capsule generation algorithm.
 * Adapted from Paul Bourke's C implementation found here:
 * http://paulbourke.net/geometry/capsule/
 */

#include <unordered_map>
#include <vector>

#include "MathDefs.h"
#include "SolidMesh.h"

class Capsule : public SolidMesh {
 public:
  Capsule(int N, const scalar& radius, const scalar& halfheight);
};

#endif
