//
// This file is part of the libWetCloth open source project
//
// Copyright 2016 Gabriel Cirio
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __SHELL_BENDING_FORCE_H__
#define __SHELL_BENDING_FORCE_H__

#include <Eigen/Core>

#include "../../Force.h"
#include "../../MathDefs.h"

class ShellBendingForce : public Force {
 protected:
  const VectorXs& m_pos;
  const VectorXs& m_rest_pos;
  const MatrixXi& m_F;
  const VectorXs& m_triangle_rest_area;

  const MatrixXi& m_E_unique;

  const MatrixXi& m_per_unique_edge_triangles;
  const MatrixXi& m_per_unique_edge_triangles_local_corners;
  const MatrixXi& m_per_triangles_unique_edges;

  const scalar& m_young_modulus;
  const scalar& m_viscous_modulus;
  const scalar& m_poisson_ratio;

  scalar m_thickness;

  bool m_apply_viscous;

  scalar m_bending_stiffness;
  scalar m_viscous_stiffness;
  VectorXs m_per_edge_rest_phi;  // theta or tantheta, depending on choice of
                                 // bending formulation
  VectorXs m_per_edge_start_phi;
  std::vector<int> m_unique_edge_usable;
  VectorXs m_multipliers;

  void computeBendingRestPhi(const VectorXs& rest_pos, VectorXs& rest_phis);

  int m_bending_mode;

 public:
  ShellBendingForce(const VectorXs& pos, const VectorXs& rest_pos,
                    const MatrixXi& F, const VectorXs& triangle_rest_area,
                    const MatrixXi& E_unique,
                    const MatrixXi& per_unique_edge_triangles,
                    const MatrixXi& per_unique_edge_triangles_local_corners,
                    const MatrixXi& per_triangles_unique_edges,
                    const scalar& young_modulus, const scalar& viscous_modulus,
                    const scalar& poisson_ratio, const scalar& thickness,
                    int bending_mode, bool apply_viscous);

  virtual ~ShellBendingForce();

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

  virtual int flag() const;
};

#endif
