//
// This file is part of the libWetCloth open source project
//
// Copyright 2016 Gabriel Cirio
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ShellMembraneForce.h"

#include <iostream>

#include "../../ThreadUtils.h"

ShellMembraneForce::~ShellMembraneForce() {}

void ShellMembraneForce::perFaceNormals(const VectorXs& V, const MatrixXi& F,
                                        MatrixXs& N) {
  N.resize(F.rows(), 3);

  int Frows = F.rows();

  threadutils::for_each(0, Frows, [&](int i) {
    const Vector3s& v1 = V.segment<3>(F(i, 1) * 4) - V.segment<3>(F(i, 0) * 4);
    const Vector3s& v2 = V.segment<3>(F(i, 2) * 4) - V.segment<3>(F(i, 0) * 4);

    N.row(i) = (v1.cross(v2)).transpose();
    const scalar r = N.row(i).norm();

    if (r > 0.0) {
      N.row(i) /= r;
    }
  });
}

ShellMembraneForce::ShellMembraneForce(
    const VectorXs& rest_pos, const VectorXs& pos, const MatrixXi& F,
    const VectorXs& triangle_rest_area, const scalar& young_modulus,
    const scalar& viscous_modulus, const scalar& poisson_ratio,
    const scalar& thickness, bool apply_viscous, bool use_approx_jacobian,
    bool use_tournier_jacobian)
    : m_rest_pos(rest_pos),
      m_pos(pos),
      m_F(F),
      m_triangle_rest_area(triangle_rest_area),
      m_young_modulus(young_modulus),
      m_viscous_modulus(viscous_modulus),
      m_poisson_ratio(poisson_ratio),
      m_apply_viscous(apply_viscous),
      m_use_approx_jacobian(use_approx_jacobian),
      m_use_tournier_jacobian(use_tournier_jacobian),
      m_thickness(thickness) {
  perFaceNormals(m_rest_pos, m_F, m_triangle_normals);
  m_membrane_ru.resize(m_F.rows(), 3);
  m_membrane_rv.resize(m_F.rows(), 3);
  m_membrane_multiplier.resize(m_F.rows() * 3);
  m_viscous_multipler.resize(m_F.rows() * 3);
  m_membrane_multiplier.setZero();
  m_viscous_multipler.setZero();

  threadutils::for_each(0, (int)m_F.rows(), [&](int f) {
    const Vector3s& x0 = m_rest_pos.segment<3>(m_F(f, 0) * 4).transpose();
    const Vector3s& x1 = m_rest_pos.segment<3>(m_F(f, 1) * 4).transpose();
    const Vector3s& x2 = m_rest_pos.segment<3>(m_F(f, 2) * 4).transpose();

    const Vector3s normal = m_triangle_normals.row(f);

    // Define (u,v) coords by LOCALLY rotating the triangle to align it with the
    // Z axis
    Eigen::Quaternion<scalar> rot =
        Eigen::Quaternion<scalar>::FromTwoVectors(normal, Vector3s::UnitZ());

    // compute uv of each vertex
    Vector2s uvi = Vector2s::Zero();
    Vector2s uvj = (rot * (x1 - x0)).segment(0, 2);  // cm
    Vector2s uvk = (rot * (x2 - x0)).segment(0, 2);

    // Determine vertex weights for strain computation
    const scalar dinv =
        1. / (uvi(0) * (uvj(1) - uvk(1)) + uvj(0) * (uvk(1) - uvi(1)) +
              uvk(0) * (uvi(1) - uvj(1)));  // cm^{-2}
    m_membrane_ru.row(f) << dinv * (uvj(1) - uvk(1)), dinv * (uvk(1) - uvi(1)),
        dinv * (uvi(1) - uvj(1));  // cm^{-1}
    m_membrane_rv.row(f) << dinv * (uvk(0) - uvj(0)), dinv * (uvi(0) - uvk(0)),
        dinv * (uvj(0) - uvi(0));
  });

  m_membrane_material_tensor_base << 1, m_poisson_ratio, 0, m_poisson_ratio, 1,
      0, 0, 0, 0.5 * (1 - m_poisson_ratio);  // 0
  m_membrane_material_viscous_tensor =
      m_membrane_material_tensor_base *
      (m_viscous_modulus * thickness / (1 - m_poisson_ratio * m_poisson_ratio));
  m_membrane_material_tensor =
      m_membrane_material_tensor_base *
      (m_young_modulus * thickness /
       (1 - m_poisson_ratio * m_poisson_ratio));  // dyne/cm
}

void ShellMembraneForce::addEnergyToTotal(const VectorXs& x, const VectorXs& v,
                                          const VectorXs& m,
                                          const VectorXs& psi,
                                          const scalar& lambda, scalar& E) {
  std::cerr << "NOT IMPLEMENTED!" << std::endl;
}

void ShellMembraneForce::addGradEToTotal(const VectorXs& x, const VectorXs& v,
                                         const VectorXs& m, const VectorXs& psi,
                                         const scalar& lambda,
                                         VectorXs& gradE) {
  for (int f = 0; f < m_F.rows(); ++f) {
    const scalar psi_base =
        (psi(m_F(f, 0)) + psi(m_F(f, 1)) + psi(m_F(f, 2))) / 3.0;
    const scalar psi_coeff = pow(psi_base, lambda);

    const Vector3s& x0 = x.segment<3>(m_F(f, 0) * 4);
    const Vector3s& x1 = x.segment<3>(m_F(f, 1) * 4);
    const Vector3s& x2 = x.segment<3>(m_F(f, 2) * 4);

    const Vector3s U = x0 * m_membrane_ru(f, 0) + x1 * m_membrane_ru(f, 1) +
                       x2 * m_membrane_ru(f, 2);  // unitless
    const Vector3s V = x0 * m_membrane_rv(f, 0) + x1 * m_membrane_rv(f, 1) +
                       x2 * m_membrane_rv(f, 2);

    const Vector3s strain = Vector3s(0.5 * (U.dot(U) - 1), 0.5 * (V.dot(V) - 1),
                                     U.dot(V));                   // unitless
    const Vector3s stress = m_membrane_material_tensor * strain;  // dyne/cm

    for (int i = 0; i < 3; ++i) {
      gradE.segment<3>(m_F(f, i) * 4) +=
          psi_coeff * m_triangle_rest_area(f) *
          (stress(0) * (m_membrane_ru(f, i) * U) +
           stress(1) * (m_membrane_rv(f, i) * V) +
           stress(2) *
               (m_membrane_ru(f, i) * V + m_membrane_rv(f, i) * U));  // dyne
    }

    if (m_apply_viscous) {
      const Vector3s& sx0 = m_start_pos.segment<3>(m_F(f, 0) * 4);
      const Vector3s& sx1 = m_start_pos.segment<3>(m_F(f, 1) * 4);
      const Vector3s& sx2 = m_start_pos.segment<3>(m_F(f, 2) * 4);

      const Vector3s sU = sx0 * m_membrane_ru(f, 0) +
                          sx1 * m_membrane_ru(f, 1) + sx2 * m_membrane_ru(f, 2);
      const Vector3s sV = sx0 * m_membrane_rv(f, 0) +
                          sx1 * m_membrane_rv(f, 1) + sx2 * m_membrane_rv(f, 2);

      const Vector3s viscous_strain =
          Vector3s(0.5 * (U.dot(U) - sU.dot(sU)), 0.5 * (V.dot(V) - sV.dot(sV)),
                   U.dot(V) - sU.dot(sV));
      const Vector3s viscous_stress =
          m_membrane_material_viscous_tensor * viscous_strain;

      for (int i = 0; i < 3; ++i) {
        gradE.segment<3>(m_F(f, i) * 4) +=
            psi_coeff * m_triangle_rest_area(f) *
            (viscous_stress(0) * (m_membrane_ru(f, i) * U) +
             viscous_stress(1) * (m_membrane_rv(f, i) * V) +
             viscous_stress(2) *
                 (m_membrane_ru(f, i) * V + m_membrane_rv(f, i) * U));
      }
    }
  }
}

void ShellMembraneForce::addHessXToTotal(const VectorXs& x, const VectorXs& v,
                                         const VectorXs& m, const VectorXs& psi,
                                         const scalar& lambda, TripletXs& hessE,
                                         int hessE_index, const scalar& dt) {
  threadutils::for_each(0, (int)m_F.rows(), [&](int f) {
    const scalar psi_base =
        (psi(m_F(f, 0)) + psi(m_F(f, 1)) + psi(m_F(f, 2))) / 3.0;
    const scalar psi_coeff = pow(psi_base, lambda);
    const Vector3s& x0 = x.segment<3>(m_F(f, 0) * 4);
    const Vector3s& x1 = x.segment<3>(m_F(f, 1) * 4);
    const Vector3s& x2 = x.segment<3>(m_F(f, 2) * 4);

    Vector3s U = x0 * m_membrane_ru(f, 0) + x1 * m_membrane_ru(f, 1) +
                 x2 * m_membrane_ru(f, 2);
    Vector3s V = x0 * m_membrane_rv(f, 0) + x1 * m_membrane_rv(f, 1) +
                 x2 * m_membrane_rv(f, 2);

    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = 0; j < 3; j++) {
        Matrix3s dphidx_membrane =
            (psi_coeff * m_triangle_rest_area(f)) *
            (m_membrane_material_tensor(0, 0) *
                 (m_membrane_ru(f, i) * m_membrane_ru(f, j) * U *
                  U.transpose()) +
             m_membrane_material_tensor(1, 1) *
                 (m_membrane_rv(f, i) * m_membrane_rv(f, j) * V *
                  V.transpose()) +
             m_membrane_material_tensor(2, 2) *
                 ((m_membrane_rv(f, i) * m_membrane_ru(f, j) * U *
                       V.transpose() +
                   m_membrane_ru(f, i) * m_membrane_rv(f, j) * V *
                       U.transpose()) +
                  (m_membrane_rv(f, i) * m_membrane_rv(f, j) * U *
                       U.transpose() +
                   m_membrane_ru(f, i) * m_membrane_ru(f, j) * V *
                       V.transpose())) +
             m_membrane_material_tensor(0, 1) *
                 (m_membrane_ru(f, i) * m_membrane_rv(f, j) * U *
                  V.transpose()) +
             m_membrane_material_tensor(1, 0) *
                 (m_membrane_rv(f, i) * m_membrane_ru(f, j) * V *
                  U.transpose()));

        if (!m_use_approx_jacobian) {
          if (m_use_tournier_jacobian) {
            dphidx_membrane +=
                (m_membrane_multiplier(f * 3 + 0) *
                     (m_membrane_ru(f, i) * m_membrane_ru(f, j)) +
                 m_membrane_multiplier(f * 3 + 1) *
                     (m_membrane_rv(f, i) * m_membrane_rv(f, j)) +
                 m_membrane_multiplier(f * 3 + 2) *
                     (m_membrane_ru(f, i) * m_membrane_rv(f, j) +
                      m_membrane_rv(f, i) * m_membrane_ru(f, j))) *
                Matrix3s::Identity();
          } else {
            Vector3s strain =
                Vector3s(0.5 * (U.dot(U) - 1), 0.5 * (V.dot(V) - 1), U.dot(V));
            Vector3s stress = m_membrane_material_tensor * strain;
            dphidx_membrane +=
                (psi_coeff * m_triangle_rest_area(f)) *
                (stress(0) * (m_membrane_ru(f, i) * m_membrane_ru(f, j)) +
                 stress(1) * (m_membrane_rv(f, i) * m_membrane_rv(f, j)) +
                 stress(2) * (m_membrane_ru(f, i) * m_membrane_rv(f, j) +
                              m_membrane_rv(f, i) * m_membrane_ru(f, j))) *
                Matrix3s::Identity();
          }
        }

        if (m_apply_viscous) {
          dphidx_membrane += (psi_coeff * m_triangle_rest_area(f)) *
                             (m_membrane_material_viscous_tensor(0, 0) *
                                  (m_membrane_ru(f, i) * m_membrane_ru(f, j) *
                                   U * U.transpose()) +
                              m_membrane_material_viscous_tensor(1, 1) *
                                  (m_membrane_rv(f, i) * m_membrane_rv(f, j) *
                                   V * V.transpose()) +
                              m_membrane_material_viscous_tensor(2, 2) *
                                  ((m_membrane_rv(f, i) * m_membrane_ru(f, j) *
                                        U * V.transpose() +
                                    m_membrane_ru(f, i) * m_membrane_rv(f, j) *
                                        V * U.transpose()) +
                                   (m_membrane_rv(f, i) * m_membrane_rv(f, j) *
                                        U * U.transpose() +
                                    m_membrane_ru(f, i) * m_membrane_ru(f, j) *
                                        V * V.transpose())) +
                              m_membrane_material_viscous_tensor(0, 1) *
                                  (m_membrane_ru(f, i) * m_membrane_rv(f, j) *
                                   U * V.transpose()) +
                              m_membrane_material_viscous_tensor(1, 0) *
                                  (m_membrane_rv(f, i) * m_membrane_ru(f, j) *
                                   V * U.transpose()));
          
          if (!m_use_approx_jacobian) {
            if (m_use_tournier_jacobian) {
              dphidx_membrane +=
                  (m_viscous_multipler(f * 3 + 0) *
                       (m_membrane_ru(f, i) * m_membrane_ru(f, j)) +
                   m_viscous_multipler(f * 3 + 1) *
                       (m_membrane_rv(f, i) * m_membrane_rv(f, j)) +
                   m_viscous_multipler(f * 3 + 2) *
                       (m_membrane_ru(f, i) * m_membrane_rv(f, j) +
                        m_membrane_rv(f, i) * m_membrane_ru(f, j))) *
                  Matrix3s::Identity();
            } else {
              const Vector3s& sx0 = m_start_pos.segment<3>(m_F(f, 0) * 4);
              const Vector3s& sx1 = m_start_pos.segment<3>(m_F(f, 1) * 4);
              const Vector3s& sx2 = m_start_pos.segment<3>(m_F(f, 2) * 4);

              const Vector3s sU = sx0 * m_membrane_ru(f, 0) +
                                  sx1 * m_membrane_ru(f, 1) +
                                  sx2 * m_membrane_ru(f, 2);
              const Vector3s sV = sx0 * m_membrane_rv(f, 0) +
                                  sx1 * m_membrane_rv(f, 1) +
                                  sx2 * m_membrane_rv(f, 2);

              const Vector3s viscous_strain = Vector3s(
                  0.5 * (U.dot(U) - sU.dot(sU)), 0.5 * (V.dot(V) - sV.dot(sV)),
                  U.dot(V) - sU.dot(sV));

              Vector3s viscous_stress =
                  m_membrane_material_viscous_tensor * viscous_strain;
              dphidx_membrane +=
                  (psi_coeff * m_triangle_rest_area(f)) *
                  (viscous_stress(0) *
                       (m_membrane_ru(f, i) * m_membrane_ru(f, j)) +
                   viscous_stress(1) *
                       (m_membrane_rv(f, i) * m_membrane_rv(f, j)) +
                   viscous_stress(2) *
                       (m_membrane_ru(f, i) * m_membrane_rv(f, j) +
                        m_membrane_rv(f, i) * m_membrane_ru(f, j))) *
                  Matrix3s::Identity();
            }
          }
        }

        int base_idx = f * 81 + (i * 3 + j) * 9;

        for (int s = 0; s < 3; ++s) {
          for (int r = 0; r < 3; ++r) {
            hessE[hessE_index + base_idx + s * 3 + r] = Triplets(
                m_F(f, i) * 4 + r, m_F(f, j) * 4 + s, dphidx_membrane(r, s));
          }
        }
      }
    }
  });
}

void ShellMembraneForce::updateMultipliers(
    const VectorXs& x, const VectorXs& vplus, const VectorXs& m,
    const VectorXs& psi, const scalar& lambda, const scalar& dt) {
  if (!m_use_tournier_jacobian) return;

  threadutils::for_each(0, (int)m_F.rows(), [&](int f) {
    const scalar psi_base =
        (psi(m_F(f, 0)) + psi(m_F(f, 1)) + psi(m_F(f, 2))) / 3.0;
    const scalar psi_coeff = pow(psi_base, lambda);

    const Vector3s& x0 = x.segment<3>(m_F(f, 0) * 4);
    const Vector3s& x1 = x.segment<3>(m_F(f, 1) * 4);
    const Vector3s& x2 = x.segment<3>(m_F(f, 2) * 4);

    const Vector3s& dq0 = vplus.segment<3>(m_F(f, 0) * 4);
    const Vector3s& dq1 = vplus.segment<3>(m_F(f, 1) * 4);
    const Vector3s& dq2 = vplus.segment<3>(m_F(f, 2) * 4);

    const Vector3s U = x0 * m_membrane_ru(f, 0) + x1 * m_membrane_ru(f, 1) +
                       x2 * m_membrane_ru(f, 2);  // unitless
    const Vector3s V = x0 * m_membrane_rv(f, 0) + x1 * m_membrane_rv(f, 1) +
                       x2 * m_membrane_rv(f, 2);
    const Vector3s dUdt = m_membrane_ru(f, 0) * dq0 +
                          m_membrane_ru(f, 1) * dq1 + m_membrane_ru(f, 2) * dq2;
    const Vector3s dVdt = m_membrane_rv(f, 0) * dq0 +
                          m_membrane_rv(f, 1) * dq1 + m_membrane_rv(f, 2) * dq2;

    const Vector3s strain = Vector3s(0.5 * (U.dot(U) - 1), 0.5 * (V.dot(V) - 1),
                                     U.dot(V));  // unitless

    const Vector3s Jv =
        Vector3s(U.dot(dUdt), V.dot(dVdt), U.dot(dVdt) + V.dot(dUdt));

    m_membrane_multiplier.segment<3>(f * 3) =
        m_membrane_material_tensor * psi_coeff * m_triangle_rest_area(f) *
        (strain + dt * Jv);  // dyne.cm

    if (m_apply_viscous) {
      const Vector3s& sx0 = m_start_pos.segment<3>(m_F(f, 0) * 4);
      const Vector3s& sx1 = m_start_pos.segment<3>(m_F(f, 1) * 4);
      const Vector3s& sx2 = m_start_pos.segment<3>(m_F(f, 2) * 4);

      const Vector3s sU = sx0 * m_membrane_ru(f, 0) +
                          sx1 * m_membrane_ru(f, 1) + sx2 * m_membrane_ru(f, 2);
      const Vector3s sV = sx0 * m_membrane_rv(f, 0) +
                          sx1 * m_membrane_rv(f, 1) + sx2 * m_membrane_rv(f, 2);

      const Vector3s viscous_strain =
          Vector3s(0.5 * (U.dot(U) - sU.dot(sU)), 0.5 * (V.dot(V) - sV.dot(sV)),
                   U.dot(V) - sU.dot(sV));

      m_viscous_multipler.segment<3>(f * 3) =
          m_membrane_material_viscous_tensor * psi_coeff *
          m_triangle_rest_area(f) * (viscous_strain + dt * Jv);
    }
  });
}

int ShellMembraneForce::numHessX() { return m_F.rows() * 81; }

void ShellMembraneForce::preCompute() {
  perFaceNormals(m_pos, m_F, m_triangle_normals);

  // update viscous tensor in-case dt changes
  m_membrane_material_viscous_tensor =
      m_membrane_material_tensor_base *
      (m_viscous_modulus * m_thickness /
       (1 - m_poisson_ratio * m_poisson_ratio));
}

void ShellMembraneForce::updateStartState() { m_start_pos = m_pos; }

Force* ShellMembraneForce::createNewCopy() {
  return new ShellMembraneForce(*this);
}

int ShellMembraneForce::flag() const { return 1; }
