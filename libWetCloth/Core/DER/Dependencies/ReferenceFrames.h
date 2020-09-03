//
// This file is part of the libWetCloth open source project
//
// Copyright 2012 Jean-Marie Aubry
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef REFERENCEFRAMES_H
#define REFERENCEFRAMES_H

#include "DegreesOfFreedom.h"

/**
 * \brief This class maintains the first reference frame vectors, orthogonal to
 * the edges.
 *
 * Initialized by space-parallel transportation along the rod, the reference
 * frame vectors later evolve by time-parallel transportation.
 *
 * Unit: no dimension
 */
class ReferenceFrames1 : public DependencyNode<Vec3Array> {
 public:
  ReferenceFrames1(Tangents& tangents)
      : DependencyNode<Vec3Array>(0, tangents.size()), m_tangents(tangents) {
    m_tangents.addDependent(this);

    storeInitialFrames(Vec3());
  }

  virtual const char* name() const { return "ReferenceFrames1"; }

  /**
   * \brief Propagates the initial (orthornormalized) reference frame
   * initRefFrame1 by space-parallel transportation along the rod
   * @param initRefFrame1 If zero or too far from orthogonal to the first edge,
   * is replaced by an arbitrary orthonormal vector.
   */
  void storeInitialFrames(const Vec3& initRefFrame1);

  /**
   * \brief Access to m_previousTangents.
   *
   * m_previousTangents is an internal parameter caching the tangents from the
   * previous time step so we are able to do time-parallel transportation. The
   * only need we have to access it directly is for
   * serialization/deserialization.
   */
  Vec3Array& getPreviousTangents() { return m_previousTangents; }

  void setPreviousTangents(const Vec3Array& deserialPrevTangents) {
    storeInitialFrames(Vec3());
    m_previousTangents = deserialPrevTangents;
  }

  bool checkNormality();

 protected:
  /**
   * \brief Computes new reference frames by time-parallel transportation along
   * the m_previousTangents->m_tangents motion.
   */
  virtual void compute();

  Tangents& m_tangents;
  Vec3Array m_previousTangents;
};

/**
 * \brief The second reference frame is simply maintained as tangent x
 * referenceFrames1
 */
class ReferenceFrames2 : public DependencyNode<Vec3Array> {
 public:
  ReferenceFrames2(Tangents& tangents, ReferenceFrames1& referenceFrames1)
      : DependencyNode<Vec3Array>(0, tangents.size()),
        m_tangents(tangents),
        m_referenceFrames1(referenceFrames1) {
    m_tangents.addDependent(this);
    m_referenceFrames1.addDependent(this);
  }

  virtual const char* name() const { return "ReferenceFrames2"; }

 protected:
  virtual void compute();

  Tangents& m_tangents;
  ReferenceFrames1& m_referenceFrames1;
};

/**
 * \brief This maintains the reference twist defined at the angle between the
 * (space-parallel transported) Bishop frame and the actual (time-parallel
 * transported) reference frame.
 */
class ReferenceTwists : public DependencyNode<std::vector<scalar> > {
 public:
  ReferenceTwists(Tangents& tangents, ReferenceFrames1& referenceFrames1)
      : DependencyNode<std::vector<scalar> >(1, tangents.size()),
        m_tangents(tangents),
        m_referenceFrames1(referenceFrames1) {
    m_tangents.addDependent(this);
    m_referenceFrames1.addDependent(this);
  }

  virtual const char* name() const { return "ReferenceTwists"; }

 protected:
  virtual void compute();

  Tangents& m_tangents;
  ReferenceFrames1& m_referenceFrames1;
};

#endif
