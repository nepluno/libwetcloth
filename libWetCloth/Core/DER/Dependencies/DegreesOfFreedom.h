//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DEGREESOFFREEDOM_H
#define DEGREESOFFREEDOM_H

#include "DependencyNode.h"

/**
 * Unit: cm for position dofs, no dimension for theta
 */
class DOFs : public DependencyNode<VecX> {
 public:
  DOFs(const VecX& dofValues) : DependencyNode<VecX>(dofValues) {
    assert(dofValues.size() % 4 == 3);
    m_numEdges = dofValues.size() / 4;
    setClean();
  }

  // As this class is meant to be the root of the dependency tree, we provide a
  // const getter This means that DOFs are never "computed" but updated by the
  // outer loop algorithm
  const VecX& get() const { return m_value; }

  VecX& get() { return m_value; }

  Vec3 getVertex(IndexType vtx) const {
    assert(vtx < (m_numEdges + 1));

    return get().segment<3>(4 * vtx);
  }

  void setVertex(IndexType vtx, const Vec3& point) {
    m_value.segment<3>(4 * vtx) = point;
    setDependentsDirty();
  }

  void setDof(IndexType i, const scalar& val) {
    m_value[i] = val;
    setDependentsDirty();
  }

  // Accessors to the theta degrees of freedom
  const Eigen::Map<const VecX, Eigen::Unaligned, Eigen::InnerStride<4> >
  getThetas() const {
    return Eigen::Map<const VecX, Eigen::Unaligned, Eigen::InnerStride<4> >(
        m_value.data() + 3, m_numEdges);
  }

  scalar getTheta(IndexType vtx) const {
    assert(vtx < m_numEdges);

    return get()[4 * vtx + 3];
  }

  void setThetas(const VecX& thetas, int numberOfFixedThetas = 0) {
    assert(thetas.size() == m_numEdges);

    Eigen::Map<VecX, Eigen::Unaligned, Eigen::InnerStride<4> >(
        m_value.data() + 4 * numberOfFixedThetas + 3,
        m_numEdges - numberOfFixedThetas) =
        thetas.tail(m_numEdges - numberOfFixedThetas);
    setDependentsDirty();
  }

  void setTheta(IndexType vtx, scalar theta) {
    m_value[4 * vtx + 3] = theta;
    setDependentsDirty();
  }

  IndexType getNumEdges() const { return m_numEdges; }

  IndexType getNumVertices() const { return m_numEdges + 1; }

  virtual const char* name() const { return "DOFs"; }

 protected:
  virtual void compute()  // Not implemented as this is an pure input node
  {
    std::cerr << "DegreesOfFreedom::compute() should never be called"
              << std::endl;
  }

 private:
  IndexType m_numEdges;
};

/**
 * Unit: cm
 */
class Edges : public DependencyNode<Vec3Array> {
 public:
  Edges(DOFs& dofs)
      : DependencyNode<Vec3Array>(0, dofs.getNumEdges()), m_dofs(dofs) {
    m_dofs.addDependent(this);
  }

  virtual const char* name() const { return "Edges"; }

 protected:
  virtual void compute();

  DOFs& m_dofs;
};

/**
 * Unit: cm
 */
class Lengths : public DependencyNode<std::vector<scalar> > {
 public:
  Lengths(Edges& edges)
      : DependencyNode<std::vector<scalar> >(0, edges.size()), m_edges(edges) {
    m_edges.addDependent(this);
  }

  virtual const char* name() const { return "Lengths"; }

 protected:
  virtual void compute();

  Edges& m_edges;
};

/**
 * Unit: no dimension
 */
class Tangents : public DependencyNode<Vec3Array> {
 public:
  Tangents(Edges& edges, Lengths& lengths)
      : DependencyNode<Vec3Array>(0, edges.size()),
        m_edges(edges),
        m_lengths(lengths) {
    m_edges.addDependent(this);
    m_lengths.addDependent(this);
  }

  virtual const char* name() const { return "Tangents"; }

 protected:
  virtual void compute();

  Edges& m_edges;
  Lengths& m_lengths;
};

/**
 * Unit: no dimension
 */
class CurvatureBinormals : public DependencyNode<Vec3Array> {
 public:
  CurvatureBinormals(Tangents& tangents)
      : DependencyNode<Vec3Array>(1, tangents.size()), m_tangents(tangents) {
    m_tangents.addDependent(this);
  }

  virtual const char* name() const { return "CurvatureBinormals"; }

 protected:
  virtual void compute();

  Tangents& m_tangents;
};

/**
 * Unit: no dimension
 */
class TrigThetas : public DependencyNode<std::pair<VecX, VecX> > {
 public:
  TrigThetas(DOFs& dofs)
      : DependencyNode<std::pair<VecX, VecX> >(std::make_pair(VecX(), VecX())),
        m_dofs(dofs) {
    m_dofs.addDependent(this);
  }

  virtual const char* name() const { return "TrigThetas"; }

  const VecX& getSines() { return get().first; }

  const VecX& getCosines() { return get().second; }

 protected:
  virtual void compute();

  DOFs& m_dofs;

 private:
  void vdSinCos(const int n, const double a[], double r1[], double r2[]);
};

#endif
