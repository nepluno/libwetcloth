//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "AttachForce.h"

#include "TwoDScene.h"

AttachForce::AttachForce(const int pidx,
                         const std::shared_ptr<TwoDScene>& scene,
                         const scalar& k, const scalar& k_twist,
                         const scalar& b, const scalar& b_twist)
    : Force(),
      m_pidx(pidx),
      m_scene(scene),
      m_k(k),
      m_k_twist(k_twist),
      m_b(b),
      m_b_twist(b_twist) {
  assert(pidx >= 0);
  assert(m_k >= 0.0);
}

AttachForce::~AttachForce() {}

void AttachForce::addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                   const VectorXs& m, const VectorXs& psi,
                                   const scalar& lambda, scalar& E) {
  assert(x.size() == v.size());
  assert(x.size() % 4 == 0);

  const Vector4s& endpoint = m_scene->getRestPos().segment<4>(m_pidx * 4);

  E += 0.5 * m_k *
           (x.segment<3>(4 * m_pidx) - endpoint.segment(0, 3)).squaredNorm() +
       0.5 * m_k_twist * (x(4 * m_pidx + 3) - endpoint(3)) *
           (x(4 * m_pidx + 3) - endpoint(3));

  if (m_b > 0.0 || m_b_twist > 0.0) {
    E +=
        0.5 * m_b *
            (x.segment<3>(4 * m_pidx) - m_start_x.segment(0, 3)).squaredNorm() +
        0.5 * m_b * (x(4 * m_pidx + 3) - m_start_x(3)) *
            (x(4 * m_pidx + 3) - m_start_x(3));
  }
}

int AttachForce::getParticleIndex() const { return m_pidx; }

void AttachForce::addGradEToTotal(const VectorXs& x, const VectorXs& v,
                                  const VectorXs& m, const VectorXs& psi,
                                  const scalar& lambda, VectorXs& gradE) {
  const Vector4s& endpoint = m_scene->getRestPos().segment<4>(m_pidx * 4);
  // Compute the elastic component
  gradE.segment<3>(4 * m_pidx) +=
      m_k * (x.segment<3>(4 * m_pidx) - endpoint.segment(0, 3));

  if (m_b > 0.0) {
    gradE.segment<3>(4 * m_pidx) +=
        m_b * (x.segment<3>(4 * m_pidx) - m_start_x.segment(0, 3));
  }

  gradE(4 * m_pidx + 3) += m_k_twist * (x(4 * m_pidx + 3) - endpoint(3));

  if (m_b_twist > 0.0) {
    gradE(4 * m_pidx + 3) += m_b_twist * (x(4 * m_pidx + 3) - m_start_x(3));
  }
}

void AttachForce::addHessXToTotal(const VectorXs& x, const VectorXs& v,
                                  const VectorXs& m, const VectorXs& psi,
                                  const scalar& lambda, TripletXs& hessE,
                                  int hessE_index, const scalar& dt) {
  hessE[hessE_index + 0] =
      Triplets(4 * m_pidx + 0, 4 * m_pidx + 0, (m_k + m_b));
  hessE[hessE_index + 1] =
      Triplets(4 * m_pidx + 1, 4 * m_pidx + 1, (m_k + m_b));
  hessE[hessE_index + 2] =
      Triplets(4 * m_pidx + 2, 4 * m_pidx + 2, (m_k + m_b));
}

void AttachForce::addAngularHessXToTotal(const VectorXs& x, const VectorXs& v,
                                         const VectorXs& m, const VectorXs& psi,
                                         const scalar& lambda, TripletXs& hessE,
                                         int hessE_index, const scalar& dt) {
  hessE[hessE_index] = Triplets(m_pidx, m_pidx, (m_k_twist + m_b_twist));
}

void AttachForce::updateMultipliers(const VectorXs& x, const VectorXs& vplus,
                                    const VectorXs& m, const VectorXs& psi,
                                    const scalar& lambda, const scalar& dt) {}

scalar AttachForce::getKs() const { return m_k; }

scalar AttachForce::getKt() const { return m_k_twist; }

void AttachForce::preCompute() {}

void AttachForce::updateStartState() {
  m_start_x = m_scene->getX().segment<4>(m_pidx * 4);
}

int AttachForce::numHessX() { return 3; }

int AttachForce::numAngularHessX() { return 1; }

Force* AttachForce::createNewCopy() { return new AttachForce(*this); }

int AttachForce::flag() const { return 1; }
