//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BENDINGFORCE_H
#define BENDINGFORCE_H

#include "ViscousOrNotViscous.h"

class StrandForce;
struct ElasticParameters;

template <typename ViscousT = NonViscous>
class BendingForce {
 public:
  BendingForce() {}

  virtual ~BendingForce() {}

 public:
  static const IndexType s_first =
      1;  // The first index on which this force can apply
  static const IndexType s_last = 1;  // The last index (counting from the end)

  typedef Eigen::Matrix<scalar, 11, 1> LocalForceType;  // Vec11
  typedef Eigen::Matrix<scalar, 11, 11> LocalJacobianType;
  typedef Eigen::Matrix<scalar, 4, 1> LocalMultiplierType;

  static std::string getName() { return ViscousT::getName() + "bending"; }

  static scalar localEnergy(const StrandForce& strand, const IndexType vtx);

  static void computeLocal(LocalMultiplierType& localL,
                           const StrandForce& strand, const IndexType vtx,
                           const scalar& dt);

  static void computeLocal(LocalForceType& localF, const StrandForce& strand,
                           const IndexType vtx);

  static void computeLocal(LocalJacobianType& localJ, const StrandForce& strand,
                           const IndexType vtx);

  static void addInPosition(VecX& globalForce, const IndexType vtx,
                            const LocalForceType& localForce);

  static void addInPosition(VecX& globalMultiplier, const IndexType vtx,
                            const LocalMultiplierType& localL);

  static void accumulateCurrentE(scalar& energy, StrandForce& strand);
  static void accumulateCurrentF(VecX& force, StrandForce& strand);
};

#endif
