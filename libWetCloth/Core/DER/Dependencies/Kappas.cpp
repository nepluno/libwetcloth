//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "Kappas.h"

#include "ElasticStrandUtils.h"

void Kappas::compute() {
  m_value.resize(m_size);
  const Vec3Array& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3Array& materialFrames1 = m_materialFrames1.get();
  const Vec3Array& materialFrames2 = m_materialFrames2.get();

  for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx) {
    m_value[vtx].setZero();
  }

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    const Vec3& kb = curvatureBinormals[vtx];
    const Vec3& m1e = materialFrames1[vtx - 1];
    const Vec3& m2e = materialFrames2[vtx - 1];
    const Vec3& m1f = materialFrames1[vtx];
    const Vec3& m2f = materialFrames2[vtx];

    m_value[vtx] = Vec2(0.5 * kb.dot(m2e + m2f), -0.5 * kb.dot(m1e + m1f));
  }

  setDependentsDirty();
}

void GradKappas::compute() {
  m_value.resize(m_size);
  const std::vector<scalar>& lengths = m_lengths.get();
  const Vec3Array& tangents = m_tangents.get();
  const Vec3Array& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3Array& materialFrames1 = m_materialFrames1.get();
  const Vec3Array& materialFrames2 = m_materialFrames2.get();
  const Vec2Array& kappas = m_kappas.get();

  for (IndexType vtx = 0; vtx < m_firstValidIndex; ++vtx) {
    m_value[vtx].setZero();
  }

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    GradKType& gradKappa = m_value[vtx];
    gradKappa.setZero();

    const scalar norm_e = lengths[vtx - 1];
    const scalar norm_f = lengths[vtx];

    const Vec3& te = tangents[vtx - 1];
    const Vec3& tf = tangents[vtx];

    const Vec3& m1e = materialFrames1[vtx - 1];
    const Vec3& m2e = materialFrames2[vtx - 1];
    const Vec3& m1f = materialFrames1[vtx];
    const Vec3& m2f = materialFrames2[vtx];

    scalar chi = 1.0 + te.dot(tf);

    //    assert( chi>0 );
    if (chi <= 0) {
      std::cerr << "GradKappas::compute(): "
                << " chi = " << chi << " te = " << te << " tf = " << tf
                << std::endl;
      chi = 1e-12;
    }

    const Vec3& tilde_d1 = (m1e + m1f) / chi;
    const Vec3& tilde_d2 = (m2e + m2f) / chi;

    const Vec2& kappa = kappas[vtx];

#ifdef USE_APPROX_GRAD_KAPPA
    const scalar beta0 =
        fabs(kappa[0]) < 1e-16 ? 1.0 : (kappa[0] / 2.0 / sin(kappa[0] / 2.0));
    const scalar beta1 =
        fabs(kappa[1]) < 1e-16 ? 1.0 : (kappa[1] / 2.0 / sin(kappa[1] / 2.0));

    const Vec3& Dkappa0De = 1.0 / norm_e * (tf.cross(tilde_d2)) * beta0;
    const Vec3& Dkappa0Df = 1.0 / norm_f * (-te.cross(tilde_d2)) * beta0;
    const Vec3& Dkappa1De = 1.0 / norm_e * (-tf.cross(tilde_d1)) * beta1;
    const Vec3& Dkappa1Df = 1.0 / norm_f * (te.cross(tilde_d1)) * beta1;
#else
    const Vec3& tilde_t = (te + tf) / chi;

    const Vec3& Dkappa0De =
        1.0 / norm_e * (-kappa[0] * tilde_t + tf.cross(tilde_d2));
    const Vec3& Dkappa0Df =
        1.0 / norm_f * (-kappa[0] * tilde_t - te.cross(tilde_d2));
    const Vec3& Dkappa1De =
        1.0 / norm_e * (-kappa[1] * tilde_t - tf.cross(tilde_d1));
    const Vec3& Dkappa1Df =
        1.0 / norm_f * (-kappa[1] * tilde_t + te.cross(tilde_d1));
#endif

    gradKappa.block<3, 1>(0, 0) = -Dkappa0De;
    gradKappa.block<3, 1>(4, 0) = Dkappa0De - Dkappa0Df;
    gradKappa.block<3, 1>(8, 0) = Dkappa0Df;
    gradKappa.block<3, 1>(0, 1) = -Dkappa1De;
    gradKappa.block<3, 1>(4, 1) = Dkappa1De - Dkappa1Df;
    gradKappa.block<3, 1>(8, 1) = Dkappa1Df;

    const Vec3& kb = curvatureBinormals[vtx];

    gradKappa(3, 0) = -0.5 * kb.dot(m1e);
    gradKappa(7, 0) = -0.5 * kb.dot(m1f);
    gradKappa(3, 1) = -0.5 * kb.dot(m2e);
    gradKappa(7, 1) = -0.5 * kb.dot(m2f);
  }

  setDependentsDirty();
}

void HessKappas::compute() {
  m_value.resize(m_size);
  const std::vector<scalar>& lengths = m_lengths.get();
  const Vec3Array& tangents = m_tangents.get();
  const Vec3Array& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3Array& materialFrames1 = m_materialFrames1.get();
  const Vec3Array& materialFrames2 = m_materialFrames2.get();
  const Vec2Array& kappas = m_kappas.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    Mat11Pair& HessKappa = m_value[vtx];

    Mat11& DDkappa1 = HessKappa.first;
    Mat11& DDkappa2 = HessKappa.second;

    DDkappa1.setZero();
    DDkappa2.setZero();

    const scalar norm_e = lengths[vtx - 1];
    const scalar norm_f = lengths[vtx];
    const scalar norm2_e = square(
        norm_e);  // That's bloody stupid, taking the square of a square root.
    const scalar norm2_f = square(norm_f);

    const Vec3& te = tangents[vtx - 1];
    const Vec3& tf = tangents[vtx];

    const Vec3& m1e = materialFrames1[vtx - 1];
    const Vec3& m2e = materialFrames2[vtx - 1];
    const Vec3& m1f = materialFrames1[vtx];
    const Vec3& m2f = materialFrames2[vtx];

    scalar chi = 1.0 + te.dot(tf);

    //    assert( chi>0 );
    if (chi <= 0) {
      std::cerr << "HessKappas::compute(): "
                << " chi = " << chi << " te = " << te << " tf = " << tf
                << std::endl;
      chi = 1e-12;
    }

    const Vec3& tilde_t = (te + tf) / chi;
    const Vec3& tilde_d1 = (m1e + m1f) / chi;
    const Vec3& tilde_d2 = (m2e + m2f) / chi;

    const Vec2& kappa = kappas[vtx];

    const Vec3& kb = curvatureBinormals[vtx];

    const Mat3& tt_o_tt = outerProd<3>(tilde_t, tilde_t);
    const Mat3& tf_c_d2t_o_tt = outerProd<3>(tf.cross(tilde_d2), tilde_t);
    const Mat3& tt_o_tf_c_d2t = tf_c_d2t_o_tt.transpose();
    const Mat3& kb_o_d2e = outerProd<3>(kb, m2e);
    const Mat3& d2e_o_kb = kb_o_d2e.transpose();

    const Mat3& Id = Mat3::Identity();

    {
      const Mat3& D2kappa1De2 =
          1.0 / norm2_e *
              (2 * kappa[0] * tt_o_tt - (tf_c_d2t_o_tt + tt_o_tf_c_d2t)) -
          kappa[0] / (chi * norm2_e) * (Id - outerProd<3>(te, te)) +
          1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

      const Mat3& te_c_d2t_o_tt = outerProd<3>(te.cross(tilde_d2), tilde_t);
      const Mat3& tt_o_te_c_d2t = te_c_d2t_o_tt.transpose();
      const Mat3& kb_o_d2f = outerProd<3>(kb, m2f);
      const Mat3& d2f_o_kb = kb_o_d2f.transpose();

      const Mat3& D2kappa1Df2 =
          1.0 / norm2_f *
              (2 * kappa[0] * tt_o_tt + (te_c_d2t_o_tt + tt_o_te_c_d2t)) -
          kappa[0] / (chi * norm2_f) * (Id - outerProd<3>(tf, tf)) +
          1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

      const Mat3& D2kappa1DeDf =
          -kappa[0] / (chi * norm_e * norm_f) * (Id + outerProd<3>(te, tf)) +
          1.0 / (norm_e * norm_f) *
              (2 * kappa[0] * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t -
               crossMat(tilde_d2));
      const Mat3& D2kappa1DfDe = D2kappa1DeDf.transpose();

      const scalar D2kappa1Dthetae2 = -0.5 * kb.dot(m2e);
      const scalar D2kappa1Dthetaf2 = -0.5 * kb.dot(m2f);
      const Vec3& D2kappa1DeDthetae =
          1.0 / norm_e *
          (0.5 * kb.dot(m1e) * tilde_t - 1.0 / chi * tf.cross(m1e));
      const Vec3& D2kappa1DeDthetaf =
          1.0 / norm_e *
          (0.5 * kb.dot(m1f) * tilde_t - 1.0 / chi * tf.cross(m1f));
      const Vec3& D2kappa1DfDthetae =
          1.0 / norm_f *
          (0.5 * kb.dot(m1e) * tilde_t + 1.0 / chi * te.cross(m1e));
      const Vec3& D2kappa1DfDthetaf =
          1.0 / norm_f *
          (0.5 * kb.dot(m1f) * tilde_t + 1.0 / chi * te.cross(m1f));

      DDkappa1.block<3, 3>(0, 0) = D2kappa1De2;
      DDkappa1.block<3, 3>(0, 4) = -D2kappa1De2 + D2kappa1DeDf;
      DDkappa1.block<3, 3>(4, 0) = -D2kappa1De2 + D2kappa1DfDe;
      DDkappa1.block<3, 3>(4, 4) =
          D2kappa1De2 - (D2kappa1DeDf + D2kappa1DfDe) + D2kappa1Df2;
      DDkappa1.block<3, 3>(0, 8) = -D2kappa1DeDf;
      DDkappa1.block<3, 3>(8, 0) = -D2kappa1DfDe;
      DDkappa1.block<3, 3>(4, 8) = D2kappa1DeDf - D2kappa1Df2;
      DDkappa1.block<3, 3>(8, 4) = D2kappa1DfDe - D2kappa1Df2;
      DDkappa1.block<3, 3>(8, 8) = D2kappa1Df2;
      DDkappa1(3, 3) = D2kappa1Dthetae2;
      DDkappa1(7, 7) = D2kappa1Dthetaf2;
      DDkappa1(3, 7) = DDkappa1(7, 3) = 0.;
      DDkappa1.block<3, 1>(0, 3) = -D2kappa1DeDthetae;
      DDkappa1.block<1, 3>(3, 0) = DDkappa1.block<3, 1>(0, 3).transpose();
      DDkappa1.block<3, 1>(4, 3) = D2kappa1DeDthetae - D2kappa1DfDthetae;
      DDkappa1.block<1, 3>(3, 4) = DDkappa1.block<3, 1>(4, 3).transpose();
      DDkappa1.block<3, 1>(8, 3) = D2kappa1DfDthetae;
      DDkappa1.block<1, 3>(3, 8) = DDkappa1.block<3, 1>(8, 3).transpose();
      DDkappa1.block<3, 1>(0, 7) = -D2kappa1DeDthetaf;
      DDkappa1.block<1, 3>(7, 0) = DDkappa1.block<3, 1>(0, 7).transpose();
      DDkappa1.block<3, 1>(4, 7) = D2kappa1DeDthetaf - D2kappa1DfDthetaf;
      DDkappa1.block<1, 3>(7, 4) = DDkappa1.block<3, 1>(4, 7).transpose();
      DDkappa1.block<3, 1>(8, 7) = D2kappa1DfDthetaf;
      DDkappa1.block<1, 3>(7, 8) = DDkappa1.block<3, 1>(8, 7).transpose();

      assert(isSymmetric(DDkappa1));
    }

    {
      const Mat3& tf_c_d1t_o_tt = outerProd<3>(tf.cross(tilde_d1), tilde_t);
      const Mat3& tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
      const Mat3& kb_o_d1e = outerProd<3>(kb, m1e);
      const Mat3& d1e_o_kb = kb_o_d1e.transpose();

      const Mat3& D2kappa2De2 =
          1.0 / norm2_e *
              (2 * kappa[1] * tt_o_tt + (tf_c_d1t_o_tt + tt_o_tf_c_d1t)) -
          kappa[1] / (chi * norm2_e) * (Id - outerProd<3>(te, te)) -
          1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

      const Mat3& te_c_d1t_o_tt = outerProd<3>(te.cross(tilde_d1), tilde_t);
      const Mat3& tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
      const Mat3& kb_o_d1f = outerProd<3>(kb, m1f);
      const Mat3& d1f_o_kb = kb_o_d1f.transpose();

      const Mat3& D2kappa2Df2 =
          1.0 / norm2_f *
              (2 * kappa[1] * tt_o_tt - (te_c_d1t_o_tt + tt_o_te_c_d1t)) -
          kappa[1] / (chi * norm2_f) * (Id - outerProd<3>(tf, tf)) -
          1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

      const Mat3& D2kappa2DeDf =
          -kappa[1] / (chi * norm_e * norm_f) * (Id + outerProd<3>(te, tf)) +
          1.0 / (norm_e * norm_f) *
              (2 * kappa[1] * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t +
               crossMat(tilde_d1));
      const Mat3& D2kappa2DfDe = D2kappa2DeDf.transpose();

      const scalar D2kappa2Dthetae2 = 0.5 * kb.dot(m1e);
      const scalar D2kappa2Dthetaf2 = 0.5 * kb.dot(m1f);
      const Vec3& D2kappa2DeDthetae =
          1.0 / norm_e *
          (0.5 * kb.dot(m2e) * tilde_t - 1.0 / chi * tf.cross(m2e));
      const Vec3& D2kappa2DeDthetaf =
          1.0 / norm_e *
          (0.5 * kb.dot(m2f) * tilde_t - 1.0 / chi * tf.cross(m2f));
      const Vec3& D2kappa2DfDthetae =
          1.0 / norm_f *
          (0.5 * kb.dot(m2e) * tilde_t + 1.0 / chi * te.cross(m2e));
      const Vec3& D2kappa2DfDthetaf =
          1.0 / norm_f *
          (0.5 * kb.dot(m2f) * tilde_t + 1.0 / chi * te.cross(m2f));

      DDkappa2.block<3, 3>(0, 0) = D2kappa2De2;
      DDkappa2.block<3, 3>(0, 4) = -D2kappa2De2 + D2kappa2DeDf;
      DDkappa2.block<3, 3>(4, 0) = -D2kappa2De2 + D2kappa2DfDe;
      DDkappa2.block<3, 3>(4, 4) =
          D2kappa2De2 - (D2kappa2DeDf + D2kappa2DfDe) + D2kappa2Df2;
      DDkappa2.block<3, 3>(0, 8) = -D2kappa2DeDf;
      DDkappa2.block<3, 3>(8, 0) = -D2kappa2DfDe;
      DDkappa2.block<3, 3>(4, 8) = D2kappa2DeDf - D2kappa2Df2;
      DDkappa2.block<3, 3>(8, 4) = D2kappa2DfDe - D2kappa2Df2;
      DDkappa2.block<3, 3>(8, 8) = D2kappa2Df2;
      DDkappa2(3, 3) = D2kappa2Dthetae2;
      DDkappa2(7, 7) = D2kappa2Dthetaf2;
      DDkappa2(3, 7) = DDkappa2(7, 3) = 0.;
      DDkappa2.block<3, 1>(0, 3) = -D2kappa2DeDthetae;
      DDkappa2.block<1, 3>(3, 0) = DDkappa2.block<3, 1>(0, 3).transpose();
      DDkappa2.block<3, 1>(4, 3) = D2kappa2DeDthetae - D2kappa2DfDthetae;
      DDkappa2.block<1, 3>(3, 4) = DDkappa2.block<3, 1>(4, 3).transpose();
      DDkappa2.block<3, 1>(8, 3) = D2kappa2DfDthetae;
      DDkappa2.block<1, 3>(3, 8) = DDkappa2.block<3, 1>(8, 3).transpose();
      DDkappa2.block<3, 1>(0, 7) = -D2kappa2DeDthetaf;
      DDkappa2.block<1, 3>(7, 0) = DDkappa2.block<3, 1>(0, 7).transpose();
      DDkappa2.block<3, 1>(4, 7) = D2kappa2DeDthetaf - D2kappa2DfDthetaf;
      DDkappa2.block<1, 3>(7, 4) = DDkappa2.block<3, 1>(4, 7).transpose();
      DDkappa2.block<3, 1>(8, 7) = D2kappa2DfDthetaf;
      DDkappa2.block<1, 3>(7, 8) = DDkappa2.block<3, 1>(8, 7).transpose();

      assert(isSymmetric(DDkappa2));
    }
  }

  setDependentsDirty();
}

void ThetaHessKappas::compute() {
  m_value.resize(m_size);
  const Vec3Array& curvatureBinormals = m_curvatureBinormals.get();
  const Vec3Array& materialFrames1 = m_materialFrames1.get();
  const Vec3Array& materialFrames2 = m_materialFrames2.get();

  for (IndexType vtx = m_firstValidIndex; vtx < size(); ++vtx) {
    ThetaHessKType& HessKappa = m_value[vtx];
    Mat2& DDkappa1 = HessKappa.first;
    Mat2& DDkappa2 = HessKappa.second;

    DDkappa1.setZero();
    DDkappa2.setZero();

    const Vec3& m1e = materialFrames1[vtx - 1];
    const Vec3& m2e = materialFrames2[vtx - 1];
    const Vec3& m1f = materialFrames1[vtx];
    const Vec3& m2f = materialFrames2[vtx];
    const Vec3& kb = curvatureBinormals[vtx];

    const scalar D2kappa1Dthetae2 = -0.5 * kb.dot(m2e);
    const scalar D2kappa1Dthetaf2 = -0.5 * kb.dot(m2f);
    const scalar D2kappa2Dthetae2 = 0.5 * kb.dot(m1e);
    const scalar D2kappa2Dthetaf2 = 0.5 * kb.dot(m1f);

    DDkappa1(0, 0) = D2kappa1Dthetae2;
    DDkappa1(1, 1) = D2kappa1Dthetaf2;
    DDkappa1(0, 1) = DDkappa1(1, 0) = 0;

    DDkappa2(0, 0) = D2kappa2Dthetae2;
    DDkappa2(1, 1) = D2kappa2Dthetaf2;
    DDkappa2(0, 1) = DDkappa2(1, 0) = 0;
  }

  setDependentsDirty();
}
