//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifdef WIN32
#define NOMINMAX
#endif

#include "DistanceFields.h"

#include <numeric>

#include "Capsule.h"
#include "Icosphere.h"
#include "MathUtilities.h"
#include "RoundCornerBox.h"
#include "RoundCylinder.h"
#include "TwoDScene.h"
#include "array3.h"
#include "array3_utils.h"
#include "makelevelset3.h"
#include "sorter.h"

using namespace mathutils;

inline scalar sphere_phi(const Vector3s& position, const Vector3s& centre,
                         scalar radius) {
  return ((position - centre).norm() - radius);
}

inline void sphere_phi_bbx(const Vector3s& centre, scalar radius,
                           Vector3s& bbx_low, Vector3s& bbx_high) {
  bbx_low = centre - Vector3s(radius, radius, radius);
  bbx_high = centre + Vector3s(radius, radius, radius);
}

inline scalar capsule_phi(const Vector3s& position, const Vector3s& centre,
                          const scalar& radius, const scalar& halflength) {
  Vector3s a = centre - Vector3s(halflength, 0, 0);
  Vector3s pa = position - a;
  Vector3s ba = Vector3s(2.0 * halflength, 0, 0);
  scalar h =
      mathutils::clamp(pa.dot(ba) / (4.0 * halflength * halflength), 0.0, 1.0);
  return (pa - ba * h).norm() - radius;
}

inline void capsule_phi_bbx(const Vector3s& centre,
                            const Eigen::Quaternion<scalar>& rot,
                            const scalar& radius, const scalar& halflength,
                            Vector3s& bbx_low, Vector3s& bbx_high) {
  const Vector3s p0 = rot * Vector3s(-halflength, 0, 0) + centre;
  const Vector3s p1 = rot * Vector3s(halflength, 0, 0) + centre;

  bbx_low = Vector3s(std::min(p0(0), p1(0)), std::min(p0(1), p1(1)),
                     std::min(p0(2), p1(2))) -
            Vector3s(radius, radius, radius);
  bbx_high = Vector3s(std::max(p0(0), p1(0)), std::max(p0(1), p1(1)),
                      std::max(p0(2), p1(2))) +
             Vector3s(radius, radius, radius);
}

inline scalar box_phi(const Vector3s& position, const Vector3s& centre,
                      const Vector3s& expand, const scalar& radius) {
  scalar dx = fabs(position[0] - centre[0]) - expand[0];
  scalar dy = fabs(position[1] - centre[1]) - expand[1];
  scalar dz = fabs(position[2] - centre[2]) - expand[2];
  scalar dax = std::max(dx, 0.0);
  scalar day = std::max(dy, 0.0);
  scalar daz = std::max(dz, 0.0);
  return std::min(std::max(std::max(dx, dy), dz), 0.0) +
         sqrt(dax * dax + day * day + daz * daz) - radius;
}

inline void box_phi_bbx(const Vector3s& centre,
                        const Eigen::Quaternion<scalar>& rot,
                        const Vector3s& expand, const scalar& radius,
                        Vector3s& bbx_low, Vector3s& bbx_high) {
  bbx_low = bbx_high = centre;

  for (int r = 0; r < 2; ++r)
    for (int s = 0; s < 2; ++s)
      for (int t = 0; t < 2; ++t) {
        Vector3s p = rot * Vector3s(r ? expand(0) : -expand(0),
                                    s ? expand(1) : -expand(1),
                                    t ? expand(2) : -expand(2)) +
                     centre;
        bbx_low =
            Vector3s(std::min(bbx_low(0), p(0)), std::min(bbx_low(1), p(1)),
                     std::min(bbx_low(2), p(2)));
        bbx_high =
            Vector3s(std::max(bbx_high(0), p(0)), std::max(bbx_high(1), p(1)),
                     std::max(bbx_high(2), p(2)));
      }

  bbx_low -= Vector3s(radius, radius, radius);
  bbx_high += Vector3s(radius, radius, radius);
}

inline scalar cylinder_phi(const Vector3s& position, const Vector3s& centre,
                           const scalar& radius_ext, const scalar& radius_cor,
                           const scalar& h) {
  Vector3s p = position - centre;
  Vector2s d = Vector2s(Vector2s(p(0), p(2)).norm() - radius_ext + radius_cor,
                        fabs(p(1)) - h);
  return std::min(std::max(d(0), d(1)), 0.0) +
         Vector2s(std::max(d(0), 0.0), std::max(d(1), 0.0)).norm() - radius_cor;
}

inline void cylinder_phi_bbx(const Vector3s& centre,
                             const Eigen::Quaternion<scalar>& rot,
                             const scalar& radius_ext, const scalar& radius_cor,
                             const scalar& h, Vector3s& bbx_low,
                             Vector3s& bbx_high) {
  Vector3s ptb = rot * Vector3s(0, h + radius_cor, 0) + centre;
  Vector3s pta = rot * Vector3s(0, -h - radius_cor, 0) + centre;
  Vector3s a = ptb - pta;
  scalar da = a.dot(a);

  Vector3s db =
      (radius_ext + radius_cor) * Vector3s(sqrt(1.0 - a(0) * a(0) / da),
                                           sqrt(1.0 - a(1) * a(1) / da),
                                           sqrt(1.0 - a(2) * a(2) / da));
  bbx_low = Vector3s(std::min(pta(0), ptb(0)), std::min(pta(1), ptb(1)),
                     std::min(pta(2), ptb(2))) -
            db;
  bbx_high = Vector3s(std::max(pta(0), ptb(0)), std::max(pta(1), ptb(1)),
                      std::max(pta(2), ptb(2))) +
             db;

  Vector3s expansion = (bbx_high - bbx_low) * 0.5;
  bbx_low = centre - expansion;
  bbx_high = centre + expansion;
}

DistanceField::DistanceField(DISTANCE_FIELD_TYPE type_,
                             DISTANCE_FIELD_USAGE usage_, int group_,
                             int params_index_, bool sampled_)
    : type(type_),
      usage(usage_),
      parent(NULL),
      group(group_),
      params_index(params_index_),
      sampled(sampled_) {}

void DistanceField::center(Vector3s& cent) const {
  Vector3s low, high;
  local_bounding_box(low, high);
  cent = (low + high) / 2.0;
}

void DistanceField::resample_internal(const std::shared_ptr<TwoDScene>& parent,
                                      const scalar& dx, const VectorXs& exist,
                                      VectorXs& additional) {
  Vector3s bbx_low = Vector3s::Zero();
  Vector3s bbx_high = Vector3s::Zero();

  local_bounding_box(bbx_low, bbx_high);

  bbx_low -= Vector3s(dx, dx, dx);
  bbx_high += Vector3s(dx, dx, dx);

  Vector3s dxyz = (bbx_high - bbx_low) / dx;
  if (dxyz(0) <= 0.0 || dxyz(1) <= 0.0 || dxyz(2) <= 0.0) return;

  Vector3i nxyz =
      Vector3i((int)ceil(dxyz(0)), (int)ceil(dxyz(1)), (int)ceil(dxyz(2)));

  const int num_init = nxyz(0) * nxyz(1) * nxyz(2);
  const int num_exist = exist.size() / 4;
  const int num_total = num_init + num_exist;

  std::vector<int> available(num_total, 0);
  VectorXs pos_init(num_total * 3);

  const int num_nodes = (nxyz(0) + 1) * (nxyz(1) + 1) * (nxyz(2) + 1);

  // init phi
  Array3s phi_init(nxyz(0) + 1, nxyz(1) + 1, nxyz(2) + 1);

  threadutils::for_each(0, num_nodes, [&](int pidx) {
    int iz = pidx / ((nxyz(0) + 1) * (nxyz(1) + 1));
    int iy = (pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1)) / (nxyz(0) + 1);
    int ix = pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1) - iy * (nxyz(0) + 1);
    Vector3s cp = bbx_low + Vector3s(ix, iy, iz) * dx;

    Vector3s vel;

    phi_init(ix, iy, iz) = compute_phi_vel(cp, vel);
  });

  // first pass of selection for new + existing particles
  threadutils::for_each(0, num_total, [&](int pidx) {
    Vector3s cp;

    if (pidx < num_init) {
      int iz = pidx / (nxyz(0) * nxyz(1));
      int iy = (pidx - iz * nxyz(0) * nxyz(1)) / nxyz(0);
      int ix = pidx - iz * nxyz(0) * nxyz(1) - iy * nxyz(0);

      cp = bbx_low + (Vector3s(ix, iy, iz) + Vector3s::Constant(0.5)) * dx +
           Vector3s(mathutils::scalarRand(-dx, dx),
                    mathutils::scalarRand(-dx, dx),
                    mathutils::scalarRand(-dx, dx));
    } else {
      int qidx = pidx - num_init;

      cp = exist.segment<3>(qidx * 4);
    }

    Vector3s vel;

    scalar dist = compute_phi_vel(cp, vel);

    if (dist > dx * sqrt(3.0)) return;

    if (parent->computePhiVel(
            cp, vel, [](const std::shared_ptr<DistanceField>& dfptr) -> bool {
              return dfptr->usage == DFU_SOLID;
            }) < 0.0)
      return;

    Vector3s grad;

    if (pidx < num_init && dist > 0.0) {
      for (int i = 0; i < 5; ++i) {
        dist = compute_phi_vel(cp, vel);

        Vector3s icp = (cp - bbx_low) / dx;

        interpolate_gradient(grad, icp, phi_init);

        if (grad.norm() > 1e-20) grad.normalize();

        cp -= grad * dist;
      }
    }

    pos_init.segment<3>(pidx * 3) = cp;
    available[pidx] = 1;
  });

  // select all particles used to compare
  std::vector<int> mapping(num_total);

  std::partial_sum(available.begin(), available.end(), mapping.begin());

  int num_parts = mapping[num_total - 1];     // all used particles
  int num_parts_new = mapping[num_init - 1];  // all newly-used particles

  VectorXs pos_selected(num_parts * 3);

  // put all used particles into buffer
  threadutils::for_each(0, num_total, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;
    pos_selected.segment<3>(idx * 3) = pos_init.segment<3>(i * 3);
  });

  available.resize(num_parts);
  for (int i = 0; i < num_parts; ++i) available[i] = 1;

  // sort all used particles in buffer
  Sorter sorter(nxyz(0), nxyz(1), nxyz(2));

  sorter.sort(num_parts, [&](int pidx, int& i, int& j, int& k) {
    i = (int)floor((pos_selected(pidx * 3 + 0) - bbx_low(0)) / dx);
    j = (int)floor((pos_selected(pidx * 3 + 1) - bbx_low(1)) / dx);
    k = (int)floor((pos_selected(pidx * 3 + 2) - bbx_low(2)) / dx);
  });

  // for each new particles, do second pass of selection by dart-picking
  sorter.for_each_bucket_particles_colored_randomized(
      [&](int pidx, int bucket_idx) {
        // we ignore existing particles
        if (pidx >= num_parts_new || !available[pidx]) return;

        const Vector3s& pos = pos_selected.segment<3>(pidx * 3);
        sorter.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int np_bucket_idx) -> bool {
              if (pidx == npidx || !available[npidx]) return false;

              const Vector3s& npos = pos_selected.segment<3>(npidx * 3);
              const scalar dist = (npos - pos).norm();

              if (dist < 0.5773502692 * dx) {
                available[pidx] = 0;
                return true;
              }

              return false;
            });
      },
      3);

  // only map new particles
  mapping.resize(num_parts_new);
  std::partial_sum(available.begin(), available.begin() + num_parts_new,
                   mapping.begin());

  int num_results = mapping[num_parts_new - 1];

  additional.resize(num_results * 3);

  // put passed new particles into output
  threadutils::for_each(0, num_parts_new, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;

    additional.segment<3>(idx * 3) = pos_selected.segment<3>(i * 3);
  });
}

void DistanceField::sample(const scalar& dx, VectorXs& result,
                           VectorXs& normals) {
  Vector3s bbx_low = Vector3s::Zero();
  Vector3s bbx_high = Vector3s::Zero();

  local_bounding_box(bbx_low, bbx_high);

  bbx_low -= Vector3s(dx, dx, dx);
  bbx_high += Vector3s(dx, dx, dx);

  Vector3s dxyz = (bbx_high - bbx_low) / dx;
  if (dxyz(0) <= 0.0 || dxyz(1) <= 0.0 || dxyz(2) <= 0.0) return;

  Vector3i nxyz =
      Vector3i((int)ceil(dxyz(0)), (int)ceil(dxyz(1)), (int)ceil(dxyz(2)));

  const int num_init = nxyz(0) * nxyz(1) * nxyz(2);
  std::vector<int> available(num_init, 0);
  VectorXs pos_init(num_init * 3);

  const int num_nodes = (nxyz(0) + 1) * (nxyz(1) + 1) * (nxyz(2) + 1);

  Array3s phi_init(nxyz(0) + 1, nxyz(1) + 1, nxyz(2) + 1);

  threadutils::for_each(0, num_nodes, [&](int pidx) {
    int iz = pidx / ((nxyz(0) + 1) * (nxyz(1) + 1));
    int iy = (pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1)) / (nxyz(0) + 1);
    int ix = pidx - iz * (nxyz(0) + 1) * (nxyz(1) + 1) - iy * (nxyz(0) + 1);
    Vector3s cp = bbx_low + Vector3s(ix, iy, iz) * dx;

    Vector3s vel;

    phi_init(ix, iy, iz) = compute_phi_vel(cp, vel);
  });

  threadutils::for_each(0, num_init, [&](int pidx) {
    int iz = pidx / (nxyz(0) * nxyz(1));
    int iy = (pidx - iz * nxyz(0) * nxyz(1)) / nxyz(0);
    int ix = pidx - iz * nxyz(0) * nxyz(1) - iy * nxyz(0);

    Vector3s cp =
        bbx_low + Vector3s(ix, iy, iz) * dx +
        Vector3s(mathutils::scalarRand(0.0, dx), mathutils::scalarRand(0.0, dx),
                 mathutils::scalarRand(0.0, dx));

    Vector3s vel;

    scalar dist = compute_phi_vel(cp, vel);

    if (fabs(dist) > dx * sqrt(3.0)) return;

    Vector3s grad;

    for (int i = 0; i < 5; ++i) {
      dist = compute_phi_vel(cp, vel);

      Vector3s icp = (cp - bbx_low) / dx;

      interpolate_gradient(grad, icp, phi_init);

      if (grad.norm() > 1e-20) grad.normalize();

      cp -= grad * dist;
    }

    pos_init.segment<3>(pidx * 3) = cp;
    available[pidx] = 1;
  });

  std::vector<int> mapping(num_init);

  std::partial_sum(available.begin(), available.end(), mapping.begin());

  int num_parts = mapping[mapping.size() - 1];

  VectorXs pos_selected(num_parts * 3);

  threadutils::for_each(0, num_init, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;
    pos_selected.segment<3>(idx * 3) = pos_init.segment<3>(i * 3);
  });

  available.resize(num_parts);
  for (int i = 0; i < num_parts; ++i) available[i] = 1;

  Sorter sorter(nxyz(0), nxyz(1), nxyz(2));

  sorter.sort(num_parts, [&](int pidx, int& i, int& j, int& k) {
    i = (int)floor((pos_selected(pidx * 3 + 0) - bbx_low(0)) / dx);
    j = (int)floor((pos_selected(pidx * 3 + 1) - bbx_low(1)) / dx);
    k = (int)floor((pos_selected(pidx * 3 + 2) - bbx_low(2)) / dx);
  });

  sorter.for_each_bucket_particles_colored_randomized(
      [&](int pidx, int bucket_idx) {
        if (!available[pidx]) return;

        const Vector3s& pos = pos_selected.segment<3>(pidx * 3);
        sorter.loop_neighbor_bucket_particles(
            bucket_idx, [&](int npidx, int np_bucket_idx) -> bool {
              if (pidx == npidx || !available[npidx]) return false;

              const Vector3s& npos = pos_selected.segment<3>(npidx * 3);
              const scalar dist = (npos - pos).norm();

              if (dist < 2.5980762116 * dx) {
                available[pidx] = 0;
                return true;
              }

              return false;
            });
      },
      3);

  mapping.resize(num_parts);
  std::partial_sum(available.begin(), available.end(), mapping.begin());

  int num_results = mapping[mapping.size() - 1];

  result.resize(num_results * 3);

  normals.resize(num_results * 3);

  threadutils::for_each(0, num_parts, [&](int i) {
    if (!available[i]) return;
    int idx = mapping[i] - 1;

    result.segment<3>(idx * 3) = pos_selected.segment<3>(i * 3);
    Vector3s icp = (result.segment<3>(idx * 3) - bbx_low) / dx;

    Vector3s grad;

    interpolate_gradient(grad, icp, phi_init);

    if (grad.norm() > 1e-20) grad.normalize();

    normals.segment<3>(idx * 3) = grad;
  });
}

bool DistanceFieldObject::check_durations(const scalar& cur_time,
                                          const scalar& cur_vol,
                                          Vector3s& shooting_vel) {
  bool ret = false;
  for (auto& dur : durations) {
    ret =
        (cur_time >= dur.start && cur_time <= dur.end) && cur_vol < dur.maxvol;
    if (ret) {
      shooting_vel = dur.vel;
      break;
    }
  }

  return ret;
}

DistanceFieldObject::DistanceFieldObject(
    const Vector3s& center_, const VectorXs& parameter_,
    DISTANCE_FIELD_TYPE type_, DISTANCE_FIELD_USAGE usage_, bool inside,
    const Vector3s& raxis, const scalar& rangle, int group_, int params_index_,
    bool sampled_, const std::vector<DF_SOURCE_DURATION>& durations_,
    const std::string& szfn, const std::string& szfn_cache)
    : DistanceField(type_, usage_, group_, params_index_, sampled_),
      center(center_),
      parameter(parameter_),
      sign(inside ? -1.0 : 1.0),
      rot(Eigen::AngleAxis<scalar>(rangle, raxis)),
      durations(durations_) {
  V.setZero();
  omega.setZero();

  switch (type_) {
    case DFT_BOX:
      mesh = std::make_shared<RoundCornerBox>(
          32, Vector3s(parameter(0), parameter(1), parameter(2)), parameter(3));
      break;
    case DFT_SPHERE:
      mesh = std::make_shared<Icosphere>(4, parameter(0));
      break;
    case DFT_CAPSULE:
      mesh = std::make_shared<Capsule>(32, parameter(0), parameter(1));
      break;
    case DFT_CYLINDER:
      mesh = std::make_shared<RoundCylinder>(32, 8, parameter(0), parameter(1),
                                             parameter(2));
      break;
    case DFT_FILE:
      mesh = std::make_shared<SolidMesh>(szfn, parameter(0));
      process_file_mesh(szfn_cache);
      break;
    default:
      mesh = nullptr;
  }

  future_center = center;
  future_rot = rot;
}

void DistanceFieldObject::resample_mesh(const scalar& dx, VectorXs& result,
                                        VectorXs& normals) {
  const int num_tris = (int)mesh->getIndices().size();
  std::vector<scalar> PDF(num_tris);

  auto& indices = mesh->getIndices();
  auto& verts = mesh->getVertices();

  scalar tot_area = 0.0;

  for (int i = 0; i < num_tris; ++i) {
    const Vector3s& v0 = verts[indices[i](0)];
    const Vector3s& v1 = verts[indices[i](1)];
    const Vector3s& v2 = verts[indices[i](2)];

    const scalar area = (v2 - v0).cross(v1 - v0).norm();

    PDF[i] = area;
    tot_area += area;
  }

  if (tot_area == 0.0) return;

  std::partial_sum(PDF.begin(), PDF.end(), PDF.begin());
  for (int i = 0; i < num_tris; ++i) {
    PDF[i] /= tot_area;
  }

  tot_area *= 0.5;

  const scalar coverage_radius =
      dx * mathutils::defaultRadiusMultiplier() * 0.5;
  const scalar coverage_area = M_PI * coverage_radius * coverage_radius;

  const int num_samples = (int)ceil(tot_area / coverage_area);

  result.resize(num_samples * 3);
  normals.resize(num_samples * 3);

  Eigen::AngleAxis<scalar> rotaa(rot);

  for (int i = 0; i < num_samples; ++i) {
    const scalar seed0 = mathutils::scalarRand(0.0, 1.0);
    std::vector<scalar>::iterator low =
        std::lower_bound(PDF.begin(), PDF.end(), seed0);
    const int tri_idx = (int)(low - PDF.begin());

    if (tri_idx < 0 || tri_idx >= num_tris) continue;

    const scalar s1 = mathutils::scalarRand(0.0, 1.0);
    const scalar s2 = mathutils::scalarRand(0.0, 1.0);

    const Vector3s& x0 = verts[indices[tri_idx](0)];
    const Vector3s& x1 = verts[indices[tri_idx](1)];
    const Vector3s& x2 = verts[indices[tri_idx](2)];

    const scalar ss1 = sqrt(s1);

    const Vector3s p =
        (1.0 - ss1) * x0 + (ss1 * (1.0 - s2)) * x1 + (s2 * ss1) * x2;

    const Vector3s n = (x2 - x0).cross(x1 - x0).normalized();

    result.segment<3>(i * 3) = rotaa * p + center;
    normals.segment<3>(i * 3) = rotaa * n * sign;
  }
}

void DistanceFieldObject::apply_global_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  future_rot = rot * future_rot;
}

void DistanceFieldObject::apply_local_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  future_rot = future_rot * rot;
}

void DistanceFieldObject::apply_translation(const Vector3s& t) {
  future_center += t;
}

scalar DistanceFieldObject::compute_phi(const Vector3s& pos) const {
  scalar phi = 0.0;

  Vector3s dx = pos - center;

  switch (type) {
    case DFT_BOX: {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * box_phi(rotp, center,
                           Vector3s(parameter(0), parameter(1), parameter(2)),
                           parameter(3));
      break;
    }

    case DFT_SPHERE: {
      phi = sign * sphere_phi(pos, center, parameter(0));
      break;
    }

    case DFT_CAPSULE: {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * capsule_phi(rotp, center, parameter(0), parameter(1));
      break;
    }

    case DFT_CYLINDER: {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp = (irot * p0 * irot.inverse()).vec() + center;
      phi = sign * cylinder_phi(rotp, center, parameter(0), parameter(1),
                                parameter(2));
      break;
    }

    case DFT_FILE: {
      Eigen::Quaternion<scalar> p0(0.0, dx(0), dx(1), dx(2));
      Eigen::Quaternion<scalar> irot = rot.conjugate();
      Vector3s rotp_coord =
          ((irot * p0 * irot.inverse()).vec() - volume_origin) / parameter(1);
      phi = sign * interpolate_value(rotp_coord, volume);
      break;
    }
    default:
      break;
  }

  return phi;
}

void DistanceFieldObject::process_file_mesh(const std::string& szfn_cache) {
  // pre-process mesh
  Vector3s centre = mesh->getCenter();
  mesh->translate(-centre);

  center += centre;

  Vector3s bbx_min, bbx_max;
  mesh->boundingBox(bbx_min, bbx_max);

  const scalar dx = parameter(1);
  bbx_min -= Vector3s::Constant(dx * 3.0);
  bbx_max += Vector3s::Constant(dx * 3.0);

  volume_origin = bbx_min;

  // check if cache exist
  if (!szfn_cache.empty()) {
    std::ifstream ifs(szfn_cache, std::ios::binary);

    if (ifs.good()) {
      // read from file directly
      read_binary_array(ifs, volume);
      ifs.close();
      return;
    }
  }

  Vector3s extend = (bbx_max - bbx_min) / dx;

  int nx = (int)ceil(extend(0));
  int ny = (int)ceil(extend(1));
  int nz = (int)ceil(extend(2));

  make_level_set3(mesh->getIndices(), mesh->getVertices(), volume_origin, dx,
                  nx, ny, nz, volume);

  if (!szfn_cache.empty()) {
    std::ofstream ofs(szfn_cache, std::ios::binary);
    write_binary_array(ofs, volume);
    ofs.close();
  }
}

scalar DistanceFieldObject::compute_phi_vel(const Vector3s& pos,
                                            Vector3s& vel) const {
  scalar phi = compute_phi(pos);

  Vector3s dx = pos - center;

  vel = V + omega.cross(dx);

  return phi;
}

void DistanceFieldObject::advance(const scalar& dt) {
  V = (future_center - center) / dt;
  Eigen::Quaternion<scalar> q = future_rot * rot.conjugate();
  scalar len = q.vec().norm();
  if (len > 0.0) {
    scalar angle = 2.0 * atan2(len, q.w());
    omega = q.vec() / len * angle / dt;
  } else {
    omega.setZero();
  }

  center = future_center;
  rot = future_rot;
}

bool DistanceFieldObject::local_bounding_box(Vector3s& bbx_low,
                                             Vector3s& bbx_high) const {
  switch (type) {
    case DFT_CAPSULE:
      capsule_phi_bbx(center, rot, parameter(0), parameter(1), bbx_low,
                      bbx_high);
      break;
    case DFT_CYLINDER:
      cylinder_phi_bbx(center, rot, parameter(0), parameter(1), parameter(2),
                       bbx_low, bbx_high);
      break;
    case DFT_SPHERE:
      sphere_phi_bbx(center, parameter(0), bbx_low, bbx_high);
      break;
    case DFT_BOX:
      box_phi_bbx(center, rot,
                  Vector3s(parameter(0), parameter(1), parameter(2)),
                  parameter(3), bbx_low, bbx_high);
      break;
    case DFT_FILE:
      mesh->boundingBox(bbx_low, bbx_high, rot);
      bbx_low += center;
      bbx_high += center;
      break;
    default:
      break;
  }

  return sign < 0.0;
}

void DistanceFieldObject::render(
    const std::function<void(const std::vector<Vector3s>&,
                             const std::vector<Vector3i>&,
                             const Eigen::Quaternion<scalar>&, const Vector3s&,
                             const scalar&)>& func) const {
  if (!mesh) return;

  func(mesh->getVertices(), mesh->getIndices(), rot, center, sign);
}

DistanceFieldOperator::DistanceFieldOperator(DISTANCE_FIELD_TYPE type_,
                                             DISTANCE_FIELD_USAGE usage_,
                                             int group_, int params_index_,
                                             bool sampled_)
    : DistanceField(type_, usage_, group_, params_index_, sampled_) {}

bool DistanceFieldOperator::check_durations(const scalar& cur_time,
                                            const scalar& cur_vol,
                                            Vector3s& shooting_vel) {
  int nb = children.size();
  for (int i = 0; i < nb; ++i) {
    if (children[i]->check_durations(cur_time, cur_vol, shooting_vel))
      return true;
  };

  return false;
}

void DistanceFieldOperator::apply_global_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  int nb = children.size();
  threadutils::for_each(
      0, nb, [&](int i) { children[i]->apply_global_rotation(rot); });
}

void DistanceFieldOperator::resample_mesh(const scalar& dx, VectorXs& result,
                                          VectorXs& normals) {
  int nb = children.size();
  threadutils::for_each(0, nb, [&](int i) {
    VectorXs tmp_result;
    VectorXs tmp_norm;
    children[i]->resample_mesh(dx, tmp_result, tmp_norm);
    const int old_size = result.size();
    result.conservativeResize(result.size() + tmp_result.size());
    normals.conservativeResize(normals.size() + tmp_norm.size());

    result.segment(old_size, tmp_result.size()) = tmp_result;
    normals.segment(old_size, tmp_norm.size()) = tmp_norm;
  });
}

void DistanceFieldOperator::apply_local_rotation(
    const Eigen::Quaternion<scalar>& rot) {
  int nb = children.size();
  threadutils::for_each(0, nb,
                        [&](int i) { children[i]->apply_local_rotation(rot); });
}

void DistanceFieldOperator::apply_translation(const Vector3s& t) {
  int nb = children.size();
  threadutils::for_each(0, nb,
                        [&](int i) { children[i]->apply_translation(t); });
}

void DistanceFieldOperator::advance(const scalar& dt) {
  int nb = children.size();
  threadutils::for_each(0, nb, [&](int i) { children[i]->advance(dt); });
}

void DistanceFieldOperator::render(
    const std::function<void(const std::vector<Vector3s>&,
                             const std::vector<Vector3i>&,
                             const Eigen::Quaternion<scalar>&, const Vector3s&,
                             const scalar&)>& func) const {
  int nb = children.size();
  for (int i = 0; i < nb; ++i) {
    children[i]->render(func);
  }
}

scalar DistanceFieldOperator::compute_phi(const Vector3s& pos) const {
  switch (DistanceField::type) {
    case DFT_UNION: {
      scalar min_phi = 1e+20;

      for (auto& child : children) {
        scalar phi = child->compute_phi(pos);
        if (phi < min_phi) {
          min_phi = phi;
        }
      }

      return min_phi;
    }
    case DFT_INTERSECT: {
      scalar max_phi = -1e+20;

      for (auto& child : children) {
        scalar phi = child->compute_phi(pos);
        if (phi > max_phi) {
          max_phi = phi;
        }
      }

      return max_phi;
    }
    default:
      return 1e+20;
  }
}

scalar DistanceFieldOperator::compute_phi_vel(const Vector3s& pos,
                                              Vector3s& vel) const {
  switch (DistanceField::type) {
    case DFT_UNION: {
      scalar min_phi = 1e+20;

      for (auto& child : children) {
        Vector3s sub_vel;
        scalar phi = child->compute_phi_vel(pos, sub_vel);
        if (phi < min_phi) {
          min_phi = phi;
          vel = sub_vel;
        }
      }

      return min_phi;
    }
    case DFT_INTERSECT: {
      scalar max_phi = -1e+20;

      for (auto& child : children) {
        Vector3s sub_vel;
        scalar phi = child->compute_phi_vel(pos, sub_vel);
        if (phi > max_phi) {
          max_phi = phi;
          vel = sub_vel;
        }
      }

      return max_phi;
    }
    default:
      vel = Vector3s::Zero();
      return 1e+20;
  }
}

int DistanceFieldOperator::vote_param_indices() {
  int m = params_index;
  int i = 0;
  for (auto& child : children) {
    int cpi = child->vote_param_indices();
    if (i == 0) {
      m = cpi;
      ++i;
    } else if (m == cpi)
      ++i;
    else
      --i;
  }

  params_index = m;

  return m;
}

DISTANCE_FIELD_USAGE DistanceFieldOperator::vote_usage() {
  DISTANCE_FIELD_USAGE m = usage;
  int i = 0;
  for (auto& child : children) {
    DISTANCE_FIELD_USAGE cpi = child->vote_usage();
    if (i == 0) {
      m = cpi;
      ++i;
    } else if (m == cpi)
      ++i;
    else
      --i;
  }

  usage = m;

  return m;
}

bool DistanceFieldOperator::vote_sampled() {
  bool m = sampled;
  int i = 0;
  for (auto& child : children) {
    bool cpi = child->vote_sampled();
    if (i == 0) {
      m = cpi;
      ++i;
    } else if (m == cpi)
      ++i;
    else
      --i;
  }

  sampled = m;

  return m;
}

bool DistanceFieldOperator::local_bounding_box(Vector3s& bbx_low,
                                               Vector3s& bbx_high) const {
  bool inside = false;

  switch (type) {
    case DFT_UNION:
      bbx_low = Vector3s::Constant(1e+20);
      bbx_high = Vector3s::Constant(-1e+20);

      for (auto& child : children) {
        Vector3s bbx_l, bbx_h;
        inside = inside || child->local_bounding_box(bbx_l, bbx_h);

        for (int r = 0; r < 3; ++r) {
          bbx_low(r) = std::min(bbx_low(r), bbx_l(r));
          bbx_high(r) = std::max(bbx_high(r), bbx_h(r));
        }
      }
      break;
    case DFT_INTERSECT:
      bbx_low = Vector3s::Constant(-1e+20);
      bbx_high = Vector3s::Constant(1e+20);

      for (auto& child : children) {
        Vector3s bbx_l, bbx_h;
        inside = inside || child->local_bounding_box(bbx_l, bbx_h);

        for (int r = 0; r < 3; ++r) {
          bbx_low(r) = std::max(bbx_low(r), bbx_l(r));
          bbx_high(r) = std::min(bbx_high(r), bbx_h(r));
        }
      }
      break;

    default:
      break;
  }

  return inside;
}
