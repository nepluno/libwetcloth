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

    m_value[vtx] = Vec4(kb.dot(m2e), -kb.dot(m1e), kb.dot(m2f), -kb.dot(m1f));
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
  const Vec4Array& kappas = m_kappas.get();

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

    const Vec3& tilde_t = (te + tf) / chi;
    const Vec3& tilde_d1e = (2.0 * m1e) / chi;
    const Vec3& tilde_d1f = (2.0 * m1f) / chi;
    const Vec3& tilde_d2e = (2.0 * m2e) / chi;
    const Vec3& tilde_d2f = (2.0 * m2f) / chi;

    const Vec4& kappa = kappas[vtx];

    const Vec3& Dkappa0eDe =
        1.0 / norm_e * (-kappa[0] * tilde_t + tf.cross(tilde_d2e));
    const Vec3& Dkappa0eDf =
        1.0 / norm_f * (-kappa[0] * tilde_t - te.cross(tilde_d2e));
    const Vec3& Dkappa1eDe =
        1.0 / norm_e * (-kappa[1] * tilde_t - tf.cross(tilde_d1e));
    const Vec3& Dkappa1eDf =
        1.0 / norm_f * (-kappa[1] * tilde_t + te.cross(tilde_d1e));

    const Vec3& Dkappa0fDe =
        1.0 / norm_e * (-kappa[2] * tilde_t + tf.cross(tilde_d2f));
    const Vec3& Dkappa0fDf =
        1.0 / norm_f * (-kappa[2] * tilde_t - te.cross(tilde_d2f));
    const Vec3& Dkappa1fDe =
        1.0 / norm_e * (-kappa[3] * tilde_t - tf.cross(tilde_d1f));
    const Vec3& Dkappa1fDf =
        1.0 / norm_f * (-kappa[3] * tilde_t + te.cross(tilde_d1f));

    gradKappa.block<3, 1>(0, 0) = -Dkappa0eDe;
    gradKappa.block<3, 1>(4, 0) = Dkappa0eDe - Dkappa0eDf;
    gradKappa.block<3, 1>(8, 0) = Dkappa0eDf;
    gradKappa.block<3, 1>(0, 1) = -Dkappa1eDe;
    gradKappa.block<3, 1>(4, 1) = Dkappa1eDe - Dkappa1eDf;
    gradKappa.block<3, 1>(8, 1) = Dkappa1eDf;

    gradKappa.block<3, 1>(0, 2) = -Dkappa0fDe;
    gradKappa.block<3, 1>(4, 2) = Dkappa0fDe - Dkappa0fDf;
    gradKappa.block<3, 1>(8, 2) = Dkappa0fDf;
    gradKappa.block<3, 1>(0, 3) = -Dkappa1fDe;
    gradKappa.block<3, 1>(4, 3) = Dkappa1fDe - Dkappa1fDf;
    gradKappa.block<3, 1>(8, 3) = Dkappa1fDf;

    const Vec3& kb = curvatureBinormals[vtx];

    gradKappa(3, 0) = -kb.dot(m1e);
    gradKappa(7, 0) = 0.0;
    gradKappa(3, 1) = -kb.dot(m2e);
    gradKappa(7, 1) = 0.0;

    gradKappa(3, 2) = 0.0;
    gradKappa(7, 2) = -kb.dot(m1f);
    gradKappa(3, 3) = 0.0;
    gradKappa(7, 3) = -kb.dot(m2f);
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
  const Vec4Array& kappas = m_kappas.get();

  for (IndexType vtx = m_firstValidIndex; vtx < m_curvatureBinormals.size(); ++vtx) {
    Mat11& DDkappa0 = m_value[vtx * 4 + 0];
    Mat11& DDkappa1 = m_value[vtx * 4 + 1];
    Mat11& DDkappa2 = m_value[vtx * 4 + 2];
    Mat11& DDkappa3 = m_value[vtx * 4 + 3];

    DDkappa0.setZero();
    DDkappa1.setZero();
    DDkappa2.setZero();
    DDkappa3.setZero();

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
    const Vec3& tilde_d1e = (2.0 * m1e) / chi;
    const Vec3& tilde_d2e = (2.0 * m2e) / chi;
    const Vec3& tilde_d1f = (2.0 * m1f) / chi;
    const Vec3& tilde_d2f = (2.0 * m2f) / chi;

    const Vec4& kappa = kappas[vtx];

    const Vec3& Dkappa0eDe =
        1.0 / norm_e * (-kappa[0] * tilde_t + tf.cross(tilde_d2e));
    const Vec3& Dkappa0eDf =
        1.0 / norm_f * (-kappa[0] * tilde_t - te.cross(tilde_d2e));
    const Vec3& Dkappa1eDe =
        1.0 / norm_e * (-kappa[1] * tilde_t - tf.cross(tilde_d1e));
    const Vec3& Dkappa1eDf =
        1.0 / norm_f * (-kappa[1] * tilde_t + te.cross(tilde_d1e));

    const Vec3& Dkappa0fDe =
        1.0 / norm_e * (-kappa[2] * tilde_t + tf.cross(tilde_d2f));
    const Vec3& Dkappa0fDf =
        1.0 / norm_f * (-kappa[2] * tilde_t - te.cross(tilde_d2f));
    const Vec3& Dkappa1fDe =
        1.0 / norm_e * (-kappa[3] * tilde_t - tf.cross(tilde_d1f));
    const Vec3& Dkappa1fDf =
        1.0 / norm_f * (-kappa[3] * tilde_t + te.cross(tilde_d1f));

    const Vec3& kb = curvatureBinormals[vtx];

    const Mat3& Id = Mat3::Identity();

    const Vec3& DchiDe = 1.0 / norm_e * (Id - outerProd<3>(te, te)) * tf;
    const Vec3& DchiDf = 1.0 / norm_f * (Id - outerProd<3>(tf, tf)) * te;

    const Mat3& DtfDf = 1.0 / norm_f * (Id - outerProd<3>(tf, tf));

    const Mat3& DttDe =
        1.0 / (chi * norm_e) *
        ((Id - outerProd<3>(te, te)) -
         outerProd<3>(tilde_t, (Id - outerProd<3>(te, te)) * tf));

    const Mat3& DttDf =
        1.0 / (chi * norm_f) *
        ((Id - outerProd<3>(tf, tf)) -
         outerProd<3>(tilde_t, (Id - outerProd<3>(tf, tf)) * te));

    // 1st Hessian
    const Mat3& D2kappa0De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa0eDe) + kappa[0] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d2e), DchiDe));

    const Mat3& D2kappa0Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa0eDf) + kappa[0] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d2e), DchiDf));

    const Mat3& D2kappa0DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa0eDf) + kappa[0] * DttDf +
         outerProd<3>(1.0 / chi * tf.cross(tilde_d2e), DchiDf) +
         crossMat(tilde_d2e) * DtfDf);

    const Mat3& D2kappa0DfDe = D2kappa0DeDf.transpose();

    const scalar D2kappa0Dthetae2 = -kb.dot(m2e);
    const scalar D2kappa0Dthetaf2 = 0.0;
    const Vec3& D2kappa0DeDthetae =
        1.0 / norm_e * (kb.dot(m1e) * tilde_t - 2.0 / chi * tf.cross(m1e));
    const Vec3& D2kappa0DeDthetaf = Vec3::Zero();
    const Vec3& D2kappa0DfDthetae =
        1.0 / norm_f * (kb.dot(m1e) * tilde_t + 2.0 / chi * te.cross(m1e));
    const Vec3& D2kappa0DfDthetaf = Vec3::Zero();

    // 2nd Hessian
    const Mat3& D2kappa1De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa1eDe) + kappa[1] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d1e), DchiDe));

    const Mat3& D2kappa1Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa1eDf) + kappa[1] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d1e), DchiDf));

    const Mat3& D2kappa1DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa1eDf) + kappa[1] * DttDf -
         outerProd<3>(1.0 / chi * tf.cross(tilde_d1e), DchiDf) -
         crossMat(tilde_d1e) * DtfDf);

    const Mat3& D2kappa1DfDe = D2kappa1DeDf.transpose();

    const scalar D2kappa1Dthetae2 = kb.dot(m1e);
    const scalar D2kappa1Dthetaf2 = 0.0;
    const Vec3& D2kappa1DeDthetae =
        1.0 / norm_e * (kb.dot(m2e) * tilde_t - 2.0 / chi * tf.cross(m2e));
    const Vec3& D2kappa1DeDthetaf = Vec3::Zero();
    const Vec3& D2kappa1DfDthetae =
        1.0 / norm_f * (kb.dot(m2e) * tilde_t + 2.0 / chi * te.cross(m2e));
    const Vec3& D2kappa1DfDthetaf = Vec3::Zero();

    // 3rd Hessian
    const Mat3& D2kappa2De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa0fDe) + kappa[2] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d2f), DchiDe));

    const Mat3& D2kappa2Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa0fDf) + kappa[2] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d2f), DchiDf));

    const Mat3& D2kappa2DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa0fDf) + kappa[2] * DttDf +
         outerProd<3>(1.0 / chi * tf.cross(tilde_d2f), DchiDf) +
         crossMat(tilde_d2f) * DtfDf);

    const Mat3& D2kappa2DfDe = D2kappa2DeDf.transpose();

    const scalar D2kappa2Dthetae2 = 0.0;
    const scalar D2kappa2Dthetaf2 = -kb.dot(m2f);
    const Vec3& D2kappa2DeDthetae = Vec3::Zero();
    const Vec3& D2kappa2DeDthetaf =
        1.0 / norm_e * (kb.dot(m1f) * tilde_t - 2.0 / chi * tf.cross(m1f));
    const Vec3& D2kappa2DfDthetae = Vec3::Zero();
    const Vec3& D2kappa2DfDthetaf =
        1.0 / norm_f * (kb.dot(m1f) * tilde_t + 2.0 / chi * te.cross(m1f));

    // 4th Hessian
    const Mat3& D2kappa3De2 =
        -1.0 / norm_e *
        symPart<3>(outerProd<3>(tilde_t + te, Dkappa1fDe) + kappa[3] * DttDe +
                   outerProd<3>(1.0 / chi * tf.cross(tilde_d1f), DchiDe));

    const Mat3& D2kappa3Df2 =
        -1.0 / norm_f *
        symPart<3>(outerProd<3>(tilde_t + tf, Dkappa1fDf) + kappa[3] * DttDf +
                   outerProd<3>(1.0 / chi * te.cross(tilde_d1f), DchiDf));

    const Mat3& D2kappa3DeDf =
        -1.0 / norm_e *
        (outerProd<3>(tilde_t, Dkappa1fDf) + kappa[3] * DttDf -
         outerProd<3>(1.0 / chi * tf.cross(tilde_d1f), DchiDf) -
         crossMat(tilde_d1f) * DtfDf);

    const Mat3& D2kappa3DfDe = D2kappa3DeDf.transpose();

    const scalar D2kappa3Dthetae2 = 0.0;
    const scalar D2kappa3Dthetaf2 = kb.dot(m1f);
    const Vec3& D2kappa3DeDthetae = Vec3::Zero();
    const Vec3& D2kappa3DeDthetaf =
        1.0 / norm_e * (kb.dot(m2f) * tilde_t - 2.0 / chi * tf.cross(m2f));
    const Vec3& D2kappa3DfDthetae = Vec3::Zero();
    const Vec3& D2kappa3DfDthetaf =
        1.0 / norm_f * (kb.dot(m2f) * tilde_t + 2.0 / chi * te.cross(m2f));

    DDkappa0.block<3, 3>(0, 0) = D2kappa0De2;
    DDkappa0.block<3, 3>(0, 4) = -D2kappa0De2 + D2kappa0DeDf;
    DDkappa0.block<3, 3>(4, 0) = -D2kappa0De2 + D2kappa0DfDe;
    DDkappa0.block<3, 3>(4, 4) =
        D2kappa0De2 - (D2kappa0DeDf + D2kappa0DfDe) + D2kappa0Df2;
    DDkappa0.block<3, 3>(0, 8) = -D2kappa0DeDf;
    DDkappa0.block<3, 3>(8, 0) = -D2kappa0DfDe;
    DDkappa0.block<3, 3>(4, 8) = D2kappa0DeDf - D2kappa0Df2;
    DDkappa0.block<3, 3>(8, 4) = D2kappa0DfDe - D2kappa0Df2;
    DDkappa0.block<3, 3>(8, 8) = D2kappa0Df2;
    DDkappa0(3, 3) = D2kappa0Dthetae2;
    DDkappa0(7, 7) = D2kappa0Dthetaf2;
    DDkappa0(3, 7) = DDkappa0(7, 3) = 0.;
    DDkappa0.block<3, 1>(0, 3) = -D2kappa0DeDthetae;
    DDkappa0.block<1, 3>(3, 0) = DDkappa0.block<3, 1>(0, 3).transpose();
    DDkappa0.block<3, 1>(4, 3) = D2kappa0DeDthetae - D2kappa0DfDthetae;
    DDkappa0.block<1, 3>(3, 4) = DDkappa0.block<3, 1>(4, 3).transpose();
    DDkappa0.block<3, 1>(8, 3) = D2kappa0DfDthetae;
    DDkappa0.block<1, 3>(3, 8) = DDkappa0.block<3, 1>(8, 3).transpose();
    DDkappa0.block<3, 1>(0, 7) = -D2kappa0DeDthetaf;
    DDkappa0.block<1, 3>(7, 0) = DDkappa0.block<3, 1>(0, 7).transpose();
    DDkappa0.block<3, 1>(4, 7) = D2kappa0DeDthetaf - D2kappa0DfDthetaf;
    DDkappa0.block<1, 3>(7, 4) = DDkappa0.block<3, 1>(4, 7).transpose();
    DDkappa0.block<3, 1>(8, 7) = D2kappa0DfDthetaf;
    DDkappa0.block<1, 3>(7, 8) = DDkappa0.block<3, 1>(8, 7).transpose();

    assert(isSymmetric(DDkappa0));

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

    DDkappa3.block<3, 3>(0, 0) = D2kappa3De2;
    DDkappa3.block<3, 3>(0, 4) = -D2kappa3De2 + D2kappa3DeDf;
    DDkappa3.block<3, 3>(4, 0) = -D2kappa3De2 + D2kappa3DfDe;
    DDkappa3.block<3, 3>(4, 4) =
        D2kappa3De2 - (D2kappa3DeDf + D2kappa3DfDe) + D2kappa3Df2;
    DDkappa3.block<3, 3>(0, 8) = -D2kappa3DeDf;
    DDkappa3.block<3, 3>(8, 0) = -D2kappa3DfDe;
    DDkappa3.block<3, 3>(4, 8) = D2kappa3DeDf - D2kappa3Df2;
    DDkappa3.block<3, 3>(8, 4) = D2kappa3DfDe - D2kappa3Df2;
    DDkappa3.block<3, 3>(8, 8) = D2kappa3Df2;
    DDkappa3(3, 3) = D2kappa3Dthetae2;
    DDkappa3(7, 7) = D2kappa3Dthetaf2;
    DDkappa3(3, 7) = DDkappa3(7, 3) = 0.;
    DDkappa3.block<3, 1>(0, 3) = -D2kappa3DeDthetae;
    DDkappa3.block<1, 3>(3, 0) = DDkappa3.block<3, 1>(0, 3).transpose();
    DDkappa3.block<3, 1>(4, 3) = D2kappa3DeDthetae - D2kappa3DfDthetae;
    DDkappa3.block<1, 3>(3, 4) = DDkappa3.block<3, 1>(4, 3).transpose();
    DDkappa3.block<3, 1>(8, 3) = D2kappa3DfDthetae;
    DDkappa3.block<1, 3>(3, 8) = DDkappa3.block<3, 1>(8, 3).transpose();
    DDkappa3.block<3, 1>(0, 7) = -D2kappa3DeDthetaf;
    DDkappa3.block<1, 3>(7, 0) = DDkappa3.block<3, 1>(0, 7).transpose();
    DDkappa3.block<3, 1>(4, 7) = D2kappa3DeDthetaf - D2kappa3DfDthetaf;
    DDkappa3.block<1, 3>(7, 4) = DDkappa3.block<3, 1>(4, 7).transpose();
    DDkappa3.block<3, 1>(8, 7) = D2kappa3DfDthetaf;
    DDkappa3.block<1, 3>(7, 8) = DDkappa3.block<3, 1>(8, 7).transpose();

    assert(isSymmetric(DDkappa3));
  }

  setDependentsDirty();
}
