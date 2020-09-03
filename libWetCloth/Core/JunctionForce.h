//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef JUNCTION_FORCE_H
#define JUNCTION_FORCE_H

#include <memory>

#include "Force.h"

class JunctionForce : public Force {
  std::vector<int> m_junctions_indices;
  std::vector<int> m_base_indices;
  std::vector<std::vector<scalar> > m_bending_coeff;
  std::vector<std::vector<int> > m_junctions_edges;
  std::vector<std::vector<scalar> > m_junction_signs;
  VectorXs m_junction_orientation;
  std::shared_ptr<TwoDScene> m_scene;

  int m_count_edges;

 public:
  JunctionForce(const std::shared_ptr<TwoDScene>&);

  virtual void addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                const VectorXs& m, const VectorXs& psi,
                                const scalar& lambda, scalar& E);

  virtual void addGradEToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, const VectorXs& psi,
                               const scalar& lambda, VectorXs& gradE);

  virtual void addHessXToTotal(const VectorXs& x, const VectorXs& v,
                               const VectorXs& m, const VectorXs& psi,
                               const scalar& lambda, TripletXs& hessE,
                               int hessE_index, const scalar& dt);

  virtual void updateMultipliers(const VectorXs& x, const VectorXs& vplus,
                                 const VectorXs& m, const VectorXs& psi,
                                 const scalar& lambda, const scalar& dt);

  virtual int numHessX();

  virtual void preCompute();

  virtual void updateStartState();

  virtual Force* createNewCopy();

  virtual bool parallelized() const;

  virtual int flag() const;
};

#endif
