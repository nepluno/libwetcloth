//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "StretchingForce.h"

#include "ViscousOrNotViscous.h"

template <typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(LocalMultiplierType& localL,
                                             const StrandForce& strand,
                                             const IndexType vtx,
                                             const scalar& dt) {
  if (strand.m_strandParams->m_useTournierJacobian) {
    const scalar ks = ViscousT::ks(strand, vtx);  // dyne
    const scalar restLength = ViscousT::ellBar(strand, vtx);

    const scalar length = strand.m_strandState->m_lengths[vtx];
    const Vec3& edge = strand.m_strandState->m_tangents[vtx];
    const scalar psi_coeff = strand.m_packing_fraction[vtx];

    localL = -ks * psi_coeff *
             (length - ViscousT::ellBar(strand, vtx) +
              dt * edge.dot(strand.m_v_plus.segment<3>((vtx + 1) * 4) -
                            strand.m_v_plus.segment<3>(vtx * 4)));  // dyne.cm  
  } else {
    localL = 0.0;
  }
}

template <typename ViscousT>
scalar StretchingForce<ViscousT>::localEnergy(const StrandForce& strand,
                                              const IndexType vtx) {
  const scalar ks = ViscousT::ks(strand, vtx);
  const scalar restLength = ViscousT::ellBar(strand, vtx);

  const scalar length = strand.m_strandState->m_lengths[vtx];
  const scalar psi_coeff = strand.m_packing_fraction[vtx];

  return 0.5 * ks * square(length / restLength - 1.0) * restLength * psi_coeff;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(
    Eigen::Matrix<scalar, 6, 1>& localF, const StrandForce& strand,
    const IndexType vtx) {
  const scalar ks = ViscousT::ks(strand, vtx);
  const scalar restLength = ViscousT::ellBar(strand, vtx);

  const scalar length = strand.m_strandState->m_lengths[vtx];
  const Vec3& edge = strand.m_strandState->m_tangents[vtx];
  const scalar psi_coeff = strand.m_packing_fraction[vtx];

  Vec3 f = ks * (length / restLength - 1.0) * edge * psi_coeff;

  localF.segment<3>(0) = f;
  localF.segment<3>(3) = -f;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::computeLocal(
    Eigen::Matrix<scalar, 6, 6>& localJ, const StrandForce& strand,
    const IndexType vtx) {
  const scalar ks = ViscousT::ks(strand, vtx);
  const scalar restLength = ViscousT::ellBar(strand, vtx);

  const scalar length = strand.m_strandState->m_lengths[vtx];
  const Vec3& edge = strand.m_strandState->m_tangents[vtx];
  const scalar psi_coeff = strand.m_packing_fraction[vtx];

  Mat3 M;
  if (strand.m_strandParams->m_useApproxJacobian) {
    M = ks *
        (std::max(0.0, 1.0 / restLength - 1.0 / length) * Mat3::Identity() +
         1.0 / restLength * (edge * edge.transpose()));  
  } else if (strand.m_strandParams->m_useTournierJacobian) {
    const scalar lambda = ViscousT::stretchingMultiplier(strand, vtx) / length;
    M = (1.0 / restLength) *
        (std::max(0.0, ks * psi_coeff + lambda) * edge * edge.transpose() -
         std::min(0.0, lambda) * Mat3::Identity());
  } else {
    M = ks * ((1.0 / restLength - 1.0 / length) * Mat3::Identity() +
              1.0 / restLength * (edge * edge.transpose()));  
  }

  localJ.block<3, 3>(0, 0) = localJ.block<3, 3>(3, 3) = -M;
  localJ.block<3, 3>(0, 3) = localJ.block<3, 3>(3, 0) = M;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::addInPosition(
    ForceVectorType& globalForce, const IndexType vtx,
    const LocalForceType& localForce) {
  globalForce.segment<3>(4 * vtx) += localForce.segment<3>(0);
  globalForce.segment<3>(4 * (vtx + 1)) += localForce.segment<3>(3);
}

template <typename ViscousT>
void StretchingForce<ViscousT>::addInPosition(
    VecX& globalMultiplier, const IndexType vtx,
    const LocalMultiplierType& localL) {
  globalMultiplier(vtx) += localL;
}

template <typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentE(scalar& energy,
                                                   StrandForce& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    energy += localEnergy(strand, vtx);
  }
}

template <typename ViscousT>
void StretchingForce<ViscousT>::accumulateCurrentF(VecX& force,
                                                   StrandForce& strand) {
  for (IndexType vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx) {
    LocalForceType localF;
    computeLocal(localF, strand, vtx);
    addInPosition(force, vtx, localF);
  }
}

template class StretchingForce<NonViscous>;
template class StretchingForce<Viscous>;
