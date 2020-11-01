//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef STRAND_FORCE_H
#define STRAND_FORCE_H

#include <Eigen/Core>
#include <iostream>
#include <unordered_set>

#include "../Force.h"
#include "../TwoDScene.h"
#include "Definitions.h"
#include "Dependencies/BendingProducts.h"
#include "Dependencies/DegreesOfFreedom.h"
#include "Dependencies/ElasticityParameters.h"
#include "Dependencies/Kappas.h"
#include "Dependencies/MaterialFrames.h"
#include "Dependencies/ReferenceFrames.h"
#include "Dependencies/Twists.h"
#include "StrandParameters.h"

struct StrandState {
  StrandState(const VecX& initDofs, BendingMatrixBase& bendingMatrixBase);

  DOFs m_dofs;
  Edges m_edges;
  Lengths m_lengths;
  Tangents m_tangents;
  mutable ReferenceFrames1 m_referenceFrames1;
  mutable ReferenceFrames2 m_referenceFrames2;
  ReferenceTwists m_referenceTwists;
  Twists m_twists;
  CurvatureBinormals m_curvatureBinormals;
  TrigThetas m_trigThetas;
  mutable MaterialFrames<1> m_materialFrames1;
  mutable MaterialFrames<2> m_materialFrames2;
  Kappas m_kappas;
  GradKappas m_gradKappas;
  GradTwists m_gradTwists;
  GradTwistsSquared m_gradTwistsSquared;
  HessKappas m_hessKappas;
  HessTwists m_hessTwists;
  BendingProducts m_bendingProducts;
};

struct StartState {  // used for Viscous updates that depend of start of step
                     // state
  StartState(const VecX& initDofs);

  DOFs m_dofs;
  Edges m_edges;
  Lengths m_lengths;
  Tangents m_tangents;
  mutable ReferenceFrames1 m_referenceFrames1;
  mutable ReferenceFrames2 m_referenceFrames2;
  ReferenceTwists m_referenceTwists;
  Twists m_twists;
  CurvatureBinormals m_curvatureBinormals;
  TrigThetas m_trigThetas;
  mutable MaterialFrames<1> m_materialFrames1;
  mutable MaterialFrames<2> m_materialFrames2;
  Kappas m_kappas;
};

class StrandForce : public Force {
 public:
  StrandForce(const std::shared_ptr<TwoDScene>& scene,
              const std::vector<int>& consecutiveVertices,
              const int& parameterIndex, int globalIndex);

  virtual ~StrandForce();

  virtual Force* createNewCopy();

  virtual void preCompute();

  virtual void updateStartState();

  virtual void updateMultipliers(const VectorXs& x, const VectorXs& vplus,
                                 const VectorXs& m, const VectorXs& psi,
                                 const scalar& lambda, const scalar& dt);

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

  virtual void addAngularHessXToTotal(const VectorXs& x, const VectorXs& v,
                                      const VectorXs& m, const VectorXs& psi,
                                      const scalar& lambda, TripletXs& hessE,
                                      int hessE_index, const scalar& dt);

  virtual int numHessX();

  virtual int numAngularHessX();

  virtual int flag() const;

  virtual const char* name();

  int getGlobalIndex() const { return m_globalIndex; }
  int getNumVertices() const { return (int)m_verts.size(); }

  Vec4Array& alterRestKappas() { return m_restKappas; }

  void updateStrandState();

  // private: //todo need to make get/set calls for twodscene in order to make
  // this private again.

  int getNumEdges() const { return (int)m_verts.size() - 1; }

  void resizeInternals();
  void freezeRestShape(unsigned begin, unsigned end, scalar damping = 0.);
  void updateRestShape(const VecX& x_restshape, scalar damping = 0.);

  void updateEverythingThatDependsOnRestLengths();

  int numConstraintNonViscous();
  int numConstraintViscous();

  void recomputeGlobal();
  void clearStored();
  void getLocalAffectedVars(int colidx,
                            std::vector<std::pair<int, int> >& vars);

  template <typename AccumulatedT>
  void accumulateQuantity(AccumulatedT& accumulated);

  void accumulateHessian(TripletXs& accumulated, TripletXs& accumulated_twist);

 private:
  //// FOSSSim related //////////////////////////////////////////////////
  std::vector<int> m_verts;  // in order root to tip
  VectorXs m_packing_fraction;

  VectorXs m_v_plus;

  VectorXs m_stretching_multipliers;
  VectorXs m_bending_multipliers;
  VectorXs m_twisting_multipliers;

  VectorXs m_viscous_stretching_multipliers;
  VectorXs m_viscous_bending_multipliers;
  VectorXs m_viscous_twisting_multipliers;

  int m_globalIndex;  // Global index amongst the hairs
  std::shared_ptr<StrandParameters> m_strandParams;
  std::shared_ptr<TwoDScene> m_scene;

  // increase memory, reduce re-computation
  scalar m_strandEnergyUpdate;
  VecX m_strandForceUpdate;
  TripletXs m_strandHessianUpdate;
  TripletXs m_strandAngularHessianUpdate;

  //// Strand State (implicitly the end of timestep state, evolved from rest
  ///config) ////////////////////////
  StrandState* m_strandState;  // future state
  StartState* m_startState;    // current state

  //// Rest shape //////////////////////////////////////////////////////
  std::vector<scalar>
      m_restLengths;  // The following four members depend on m_restLengths,
                      // which is why updateEverythingThatDependsOnRestLengths()
                      // must be called
  scalar m_totalRestLength;
  std::vector<scalar> m_VoronoiLengths;     // rest length around each vertex
  std::vector<scalar> m_invVoronoiLengths;  // their inverses
  std::vector<scalar> m_vertexMasses;
  Vec4Array m_restKappas;
  std::vector<scalar> m_restTwists;

  //// Friends /////////////////////////////////////////////////////////
  friend class Viscous;
  friend class NonViscous;
  template <typename ViscousT>
  friend class StretchingForce;
  template <typename ViscousT>
  friend class BendingForce;
  template <typename ViscousT>
  friend class TwistingForce;
};

#endif
