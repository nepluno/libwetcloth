//
// This file is part of the libWetCloth open source project
//
// Copyright 2016 Gabriel Cirio
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __SHELL_MEMBRANE_FORCE_H__
#define __SHELL_MEMBRANE_FORCE_H__

#include <Eigen/Core>

#include "../../Force.h"
#include "../../MathDefs.h"

class ShellMembraneForce : public Force {
 protected:
  const VectorXs& m_rest_pos;
  const VectorXs& m_pos;
  const MatrixXi& m_F;
  const VectorXs& m_triangle_rest_area;

  MatrixXs m_triangle_normals;

  const scalar& m_young_modulus;
  const scalar& m_viscous_modulus;
  const scalar& m_poisson_ratio;
  Matrix3s m_membrane_material_tensor_base;
  Matrix3s m_membrane_material_tensor;
  Matrix3s m_membrane_material_viscous_tensor;
  MatrixXs m_membrane_ru;
  MatrixXs m_membrane_rv;

  VectorXs m_start_pos;

  VectorXs m_membrane_multiplier;
  VectorXs m_viscous_multipler;

  bool m_apply_viscous;
  scalar m_thickness;

  void perFaceNormals(const VectorXs& V, const MatrixXi& F, MatrixXs& N);

 public:
  ShellMembraneForce(const VectorXs& rest_pos, const VectorXs& pos,
                     const MatrixXi& F, const VectorXs& triangle_rest_area,
                     const scalar& young_modulus, const scalar& viscous_modulus,
                     const scalar& poisson_ratio, const scalar& thickness,
                     bool apply_viscous);

  virtual ~ShellMembraneForce();

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
