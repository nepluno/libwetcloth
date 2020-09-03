//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TwoDSceneSerializer.h"

#include <igl/boundary_loop.h>

#include <fstream>
#include <iomanip>
#include <numeric>

#include "AttachForce.h"

const int num_cloth_edge_discretization = 6;

void serialize_pos_subprog(SerializePosPacket* packet) {
  std::ofstream ofs_pos(packet->fn_pos.c_str(), std::ios::binary);
  int bufsize = packet->m_pos.size() * sizeof(scalar);
  ofs_pos.write((const char*)&bufsize, sizeof(int));
  ofs_pos.write((const char*)(packet->m_pos.data()), bufsize);

  bufsize = packet->m_d_gauss.size() * sizeof(scalar);
  ofs_pos.write((const char*)&bufsize, sizeof(int));
  ofs_pos.write((const char*)(packet->m_d_gauss.data()), bufsize);

  ofs_pos.flush();
  ofs_pos.close();

  delete packet;
}

void serialize_subprog(SerializePacket* packet) {
  std::ofstream ofs_fluid(packet->fn_fluid.c_str());
  const int num_fp = packet->m_fluid_vertices.size();
  for (int i = 0; i < num_fp; ++i) {
    ofs_fluid << "v " << std::setprecision(8) << packet->m_fluid_vertices[i](0)
              << " " << packet->m_fluid_vertices[i](1) << " "
              << packet->m_fluid_vertices[i](2) << " "
              << packet->m_fluid_radii[i] << std::endl;
  }

  std::ofstream ofs_hair(packet->fn_hairs.c_str());
  const int num_hair_vtx = packet->m_hair_vertices.size();
  for (int i = 0; i < num_hair_vtx; ++i) {
    ofs_hair << "v " << std::setprecision(8) << packet->m_hair_vertices[i](0)
             << " " << packet->m_hair_vertices[i](1) << " "
             << packet->m_hair_vertices[i](2) << " "
             << packet->m_hair_radii[i](0) << " " << packet->m_hair_radii[i](1)
             << " " << packet->m_hair_sat[i] << " " << packet->m_hair_group[i]
             << " " << packet->m_hair_vertices_rest[i](0) << " "
             << packet->m_hair_vertices_rest[i](1) << " "
             << packet->m_hair_vertices_rest[i](2) << std::endl;
  }
  for (auto& e : packet->m_hair_indices) {
    ofs_hair << "l " << (e(0) + 1) << " " << (e(1) + 1) << std::endl;
  }

  std::ofstream ofs_cloth(packet->fn_clothes.c_str());
  const int num_cloth_vtx = packet->m_dbl_face_cloth_vertices.size();
  for (int i = 0; i < num_cloth_vtx; ++i) {
    ofs_cloth << "v " << std::setprecision(8)
              << packet->m_dbl_face_cloth_vertices[i](0) << " "
              << packet->m_dbl_face_cloth_vertices[i](1) << " "
              << packet->m_dbl_face_cloth_vertices[i](2) << " "
              << packet->m_dbl_face_cloth_sat[i] << " "
              << packet->m_dbl_face_cloth_dir[i] << " "
              << packet->m_dbl_face_cloth_group[i] << " "
              << packet->m_dbl_face_cloth_vertices_rest[i](0) << " "
              << packet->m_dbl_face_cloth_vertices_rest[i](1) << " "
              << packet->m_dbl_face_cloth_vertices_rest[i](2) << " "
              << packet->m_dbl_face_cloth_vertices_central[i](0) << " "
              << packet->m_dbl_face_cloth_vertices_central[i](1) << " "
              << packet->m_dbl_face_cloth_vertices_central[i](2) << std::endl;
  }

  for (auto& f : packet->m_dbl_face_cloth_indices) {
    ofs_cloth << "f " << (f(0) + 1) << " " << (f(1) + 1) << " " << (f(2) + 1)
              << std::endl;
  }

  std::ofstream ofs_internal(packet->fn_internal_boundaries.c_str());
  for (auto& v : packet->m_internal_vertices) {
    ofs_internal << "v " << std::setprecision(8) << v(0) << " " << v(1) << " "
                 << v(2) << std::endl;
  }
  for (auto& f : packet->m_internal_indices) {
    ofs_internal << "f " << (f(0) + 1) << " " << (f(1) + 1) << " " << (f(2) + 1)
                 << std::endl;
  }

  std::ofstream ofs_external(packet->fn_external_boundaries.c_str());
  for (auto& v : packet->m_external_vertices) {
    ofs_external << "v " << std::setprecision(8) << v(0) << " " << v(1) << " "
                 << v(2) << std::endl;
  }
  for (auto& f : packet->m_external_indices) {
    ofs_external << "f " << (f(0) + 1) << " " << (f(1) + 1) << " " << (f(2) + 1)
                 << std::endl;
  }

  std::ofstream ofs_spring(packet->fn_springs.c_str());
  for (auto& v : packet->m_attach_spring_vertices) {
    ofs_spring << "v " << std::setprecision(8) << v(0) << " " << v(1) << " "
               << v(2) << std::endl;
  }
  const int num_springs = packet->m_attach_spring_vertices.size() / 2;
  for (int i = 0; i < num_springs; ++i) {
    ofs_spring << "l " << (i * 2) << " " << (i * 2 + 1) << std::endl;
  }

  ofs_fluid.flush();
  ofs_hair.flush();
  ofs_cloth.flush();
  ofs_internal.flush();
  ofs_external.flush();
  ofs_spring.flush();

  ofs_fluid.close();
  ofs_hair.close();
  ofs_cloth.close();
  ofs_internal.close();
  ofs_external.close();
  ofs_spring.close();

  std::cout << "[Frame with " << packet->fn_fluid << " written]" << std::endl;

  delete packet;
}

void TwoDSceneSerializer::serializeScene(
    TwoDScene& scene, const std::string& fn_clothes,
    const std::string& fn_hairs, const std::string& fn_fluid,
    const std::string& fn_internal_boundaries,
    const std::string& fn_external_boundaries, const std::string& fn_springs) {
  SerializePacket* data = new SerializePacket;

  updateDoubleFaceCloth(scene, data);
  updateHairs(scene, data);
  updateFluid(scene, data);
  updateMesh(scene, data);
  updateAttachSprings(scene, data);

  data->fn_clothes = fn_clothes.c_str();
  data->fn_fluid = fn_fluid.c_str();
  data->fn_hairs = fn_hairs.c_str();
  data->fn_external_boundaries = fn_external_boundaries.c_str();
  data->fn_internal_boundaries = fn_internal_boundaries.c_str();
  data->fn_springs = fn_springs.c_str();

  std::thread t(std::bind(serialize_subprog, data));
  t.detach();
}

void TwoDSceneSerializer::serializePositionOnly(TwoDScene& scene,
                                                const std::string& fn_pos) {
  SerializePosPacket* data = new SerializePosPacket;
  data->fn_pos = fn_pos.c_str();
  data->m_pos = scene.getX();
  data->m_d_gauss = scene.getGaussd();

  std::thread t(std::bind(serialize_pos_subprog, data));
  t.detach();
}

void TwoDSceneSerializer::loadPosOnly(TwoDScene& scene,
                                      std::ifstream& inputstream) {
  VectorXs& x = scene.getX();
  int bufsize;
  inputstream.read((char*)&bufsize, sizeof(int));
  inputstream.read((char*)x.data(), bufsize);

  MatrixXs& d_gauss = scene.getGaussd();
  inputstream.read((char*)&bufsize, sizeof(int));
  inputstream.read((char*)d_gauss.data(), bufsize);
}

void TwoDSceneSerializer::updateHairs(const TwoDScene& scene,
                                      SerializePacket* data) {
  const VectorXs& x = scene.getX();
  const VectorXs& rest_x = scene.getRestPos();
  const VectorXs& r = scene.getRadius();
  const VectorXs& vol = scene.getVol();
  const VectorXs& fvol = scene.getFluidVol();
  const std::vector<int>& group = scene.getParticleGroup();

  const int num_soft_elasto = scene.getNumSoftElastoParticles();
  const MatrixXi& edges = scene.getEdges();
  if (num_soft_elasto == 0 || edges.rows() == 0) return;

  std::vector<int> hair_indicator(num_soft_elasto);
  threadutils::for_each(0, num_soft_elasto, [&](int pidx) {
    const std::vector<int>& p2e = scene.getParticleEdges(pidx);

    if (p2e.size() > 0) {
      hair_indicator[pidx] = 1;
    } else {
      hair_indicator[pidx] = 0;
    }
  });

  std::partial_sum(hair_indicator.begin(), hair_indicator.end(),
                   hair_indicator.begin());

  const int num_hair_parts = hair_indicator[hair_indicator.size() - 1];

  if (num_hair_parts == 0) return;

  data->m_hair_vertices.resize(num_hair_parts);
  data->m_hair_vertices_rest.resize(num_hair_parts);
  data->m_hair_indices.resize(edges.rows());
  data->m_hair_sat.resize(num_hair_parts);
  data->m_hair_radii.resize(num_hair_parts);
  data->m_hair_group.resize(num_hair_parts);

  threadutils::for_each(0, num_soft_elasto, [&](int pidx) {
    const std::vector<int>& p2e = scene.getParticleEdges(pidx);
    if (p2e.size() == 0) return;

    const int mapped_idx = hair_indicator[pidx] - 1;
    data->m_hair_vertices[mapped_idx] = x.segment<3>(pidx * 4);
    data->m_hair_vertices_rest[mapped_idx] = rest_x.segment<3>(pidx * 4);
    data->m_hair_radii[mapped_idx] = r.segment<2>(pidx * 2);
    data->m_hair_sat[mapped_idx] = fvol(pidx) / std::max(1e-16, vol(pidx));
    data->m_hair_group[mapped_idx] = group[pidx];
  });

  threadutils::for_each(0, (int)edges.rows(), [&](int eidx) {
    const int mapped_idx_0 = hair_indicator[edges(eidx, 0)] - 1;
    const int mapped_idx_1 = hair_indicator[edges(eidx, 1)] - 1;

    data->m_hair_indices[eidx] = Vector2i(mapped_idx_0, mapped_idx_1);
  });
}

void TwoDSceneSerializer::updateFluid(const TwoDScene& scene,
                                      SerializePacket* data) {
  data->m_fluid_vertices.resize(scene.getNumFluidParticles());
  data->m_fluid_radii.resize(scene.getNumFluidParticles());

  const std::vector<int>& indices = scene.getFluidIndices();
  const VectorXs& x = scene.getX();
  const VectorXs& r = scene.getRadius();

  threadutils::for_each(0, scene.getNumFluidParticles(), [&](int idx) {
    data->m_fluid_vertices[idx] = x.segment<3>(indices[idx] * 4);
    data->m_fluid_radii[idx] = r(indices[idx] * 2 + 0);
  });
}

void TwoDSceneSerializer::updateMesh(const TwoDScene& scene,
                                     SerializePacket* data) {
  const std::vector<std::shared_ptr<DistanceField> >& fields =
      scene.getGroupDistanceField();
  for (auto& ptr : fields) {
    if (ptr->usage != DFU_SOLID) continue;

    ptr->render([&](const std::vector<Vector3s>& vertices,
                    const std::vector<Vector3i>& indices,
                    const Eigen::Quaternion<scalar>& rot,
                    const Vector3s& center, const scalar& sign) {
      if (sign > 0.0) {
        const int vert_base_idx = (int)data->m_internal_vertices.size();

        for (size_t i = 0; i < vertices.size(); ++i) {
          Eigen::Quaternion<scalar> rotp(0.0, vertices[i](0), vertices[i](1),
                                         vertices[i](2));
          Vector3s rot_x = (rot * rotp * rot.inverse()).vec() + center;

          data->m_internal_vertices.push_back(rot_x);
        }

        for (size_t i = 0; i < indices.size(); ++i) {
          data->m_internal_indices.push_back(
              indices[i] +
              Vector3i(vert_base_idx, vert_base_idx, vert_base_idx));
        }
      } else {
        int vert_base_idx = (int)data->m_external_vertices.size();

        for (size_t i = 0; i < vertices.size(); ++i) {
          Eigen::Quaternion<scalar> rotp(0.0, vertices[i](0), vertices[i](1),
                                         vertices[i](2));
          Vector3s rot_x = (rot * rotp * rot.inverse()).vec() + center;

          data->m_external_vertices.push_back(rot_x);
        }

        for (size_t i = 0; i < indices.size(); ++i) {
          Vector3i tri_inversed =
              Vector3i(indices[i](0), indices[i](2), indices[i](1));
          data->m_external_indices.push_back(
              tri_inversed +
              Vector3i(vert_base_idx, vert_base_idx, vert_base_idx));
        }

        vert_base_idx = (int)data->m_external_vertices.size();
        for (size_t i = 0; i < vertices.size(); ++i) {
          Eigen::Quaternion<scalar> rotp(0.0, vertices[i](0), vertices[i](1),
                                         vertices[i](2));
          Vector3s rot_x = (rot * rotp * rot.inverse()).vec() * 1.01 + center;

          data->m_external_vertices.push_back(rot_x);
        }

        for (size_t i = 0; i < indices.size(); ++i) {
          data->m_external_indices.push_back(
              indices[i] +
              Vector3i(vert_base_idx, vert_base_idx, vert_base_idx));
        }
      }
    });
  }
}

void TwoDSceneSerializer::initializeFaceLoops(const TwoDScene& scene) {
  const MatrixXi& faces = scene.getFaces();
  igl::boundary_loop(faces, m_face_loops);
}

void TwoDSceneSerializer::updateAttachSprings(const TwoDScene& scene,
                                              SerializePacket* data) {
  const VectorXs& x = scene.getX();
  const VectorXs& rest_x = scene.getRestPos();
  const auto& forces = scene.getAttachForces();

  data->m_attach_spring_vertices.resize(forces.size() * 2);

  const int num_f = forces.size();
  threadutils::for_each(0, num_f, [&](int i) {
    const int pidx = forces[i]->getParticleIndex();
    data->m_attach_spring_vertices[i * 2] = rest_x.segment<3>(pidx * 4);
    data->m_attach_spring_vertices[i * 2 + 1] = x.segment<3>(pidx * 4);
  });
}

void TwoDSceneSerializer::updateDoubleFaceCloth(const TwoDScene& scene,
                                                SerializePacket* data) {
  const VectorXs& x = scene.getX();
  const VectorXs& rest_x = scene.getRestPos();
  const std::vector<int> group = scene.getParticleGroup();
  const int num_edges = scene.getNumEdges();
  const int num_soft_elasto = scene.getNumSoftElastoParticles();
  const int num_faces = scene.getNumFaces();
  const MatrixXi& faces = scene.getFaces();
  const MatrixXs& gnorms = scene.getGaussNormal();
  const VectorXs& radius = scene.getRadius();
  const VectorXs& fluid_vol = scene.getFluidVol();
  const VectorXs& vol = scene.getVol();
  const VectorXs& vol_frac = scene.getVolumeFraction();

  if (num_soft_elasto == 0) return;

  // compute face norms
  std::vector<Vector3s> face_norms(num_faces);
  threadutils::for_each(0, num_faces, [&](int fidx) {
    const Vector3iT& f = faces.row(fidx);
    Vector3s t0 = (x.segment<3>(f[1] * 4) - x.segment<3>(f[0] * 4));
    Vector3s t1 = (x.segment<3>(f[2] * 4) - x.segment<3>(f[0] * 4));

    Vector3s norm = t1.cross(t0).normalized();
    face_norms[fidx] = norm;
  });

  // compute vertex norms
  std::vector<Vector3s> vert_norms(num_soft_elasto);
  std::vector<int> cloth_indicator(num_soft_elasto);

  threadutils::for_each(0, num_soft_elasto, [&](int pidx) {
    const std::vector<std::pair<int, scalar> >& p2f =
        scene.getParticleFaces(pidx);
    if (p2f.size() > 0) {
      Vector3s n = Vector3s::Zero();
      for (auto& p : p2f) {
        n += face_norms[p.first] * p.second;
      }

      vert_norms[pidx] = n.normalized();
      cloth_indicator[pidx] = 1;
    } else {
      vert_norms[pidx] = Vector3s::Zero();
      cloth_indicator[pidx] = 0;
    }
  });

  std::partial_sum(cloth_indicator.begin(), cloth_indicator.end(),
                   cloth_indicator.begin());
  const int num_single_cloth_verts =
      cloth_indicator[cloth_indicator.size() - 1];

  if (!num_single_cloth_verts) return;

  int num_loop_verts = 0;
  const int num_clothes = m_face_loops.size();
  std::vector<int> lv_base(num_clothes);

  for (int iCloth = 0; iCloth < num_clothes; ++iCloth) {
    lv_base[iCloth] = num_loop_verts;
    num_loop_verts += m_face_loops[iCloth].size();
  }

  const int total_verts = num_single_cloth_verts * 2 +
                          num_loop_verts * (num_cloth_edge_discretization - 1);
  const int total_faces =
      faces.rows() * 2 + num_loop_verts * num_cloth_edge_discretization * 2;

  data->m_dbl_face_cloth_vertices.resize(total_verts);
  data->m_dbl_face_cloth_vertices_rest.resize(total_verts);
  data->m_dbl_face_cloth_vertices_central.resize(total_verts);
  data->m_dbl_face_cloth_sat.resize(total_verts);
  data->m_dbl_face_cloth_group.resize(total_verts);
  data->m_dbl_face_cloth_dir.resize(total_verts);
  data->m_dbl_face_cloth_indices.resize(total_faces);

  threadutils::for_each(0, num_soft_elasto, [&](int pidx) {
    const std::vector<std::pair<int, scalar> >& p2f =
        scene.getParticleFaces(pidx);
    if (p2f.size() > 0) {
      const int new_pidx = cloth_indicator[pidx] - 1;
      data->m_dbl_face_cloth_vertices[new_pidx] =
          x.segment<3>(pidx * 4) + vert_norms[pidx] * radius(pidx * 2 + 0);
      data->m_dbl_face_cloth_vertices[num_single_cloth_verts + new_pidx] =
          x.segment<3>(pidx * 4) - vert_norms[pidx] * radius(pidx * 2 + 1);
      data->m_dbl_face_cloth_vertices_rest[new_pidx] =
          rest_x.segment<3>(pidx * 4);
      data->m_dbl_face_cloth_vertices_rest[num_single_cloth_verts + new_pidx] =
          rest_x.segment<3>(pidx * 4);
      data->m_dbl_face_cloth_vertices_central[new_pidx] =
          x.segment<3>(pidx * 4);
      data->m_dbl_face_cloth_vertices_central[num_single_cloth_verts +
                                              new_pidx] =
          x.segment<3>(pidx * 4);
      data->m_dbl_face_cloth_dir[new_pidx] = 0;
      data->m_dbl_face_cloth_dir[num_single_cloth_verts + new_pidx] = 2;
      data->m_dbl_face_cloth_group[new_pidx] = group[pidx];
      data->m_dbl_face_cloth_group[num_single_cloth_verts + new_pidx] =
          group[pidx];
      data->m_dbl_face_cloth_sat[new_pidx] =
          fluid_vol(pidx) / std::max(1e-16, vol(pidx));
      data->m_dbl_face_cloth_sat[num_single_cloth_verts + new_pidx] =
          fluid_vol(pidx) / std::max(1e-16, vol(pidx));
    }
  });

  threadutils::for_each(0, (int)faces.rows(), [&](int fidx) {
    const int new_pidx_0 = cloth_indicator[faces(fidx, 0)] - 1;
    const int new_pidx_1 = cloth_indicator[faces(fidx, 1)] - 1;
    const int new_pidx_2 = cloth_indicator[faces(fidx, 2)] - 1;

    data->m_dbl_face_cloth_indices[fidx] =
        Vector3i(new_pidx_0, new_pidx_2, new_pidx_1);
    data->m_dbl_face_cloth_indices[fidx + faces.rows()] =
        Vector3i(new_pidx_0 + num_single_cloth_verts,
                 new_pidx_1 + num_single_cloth_verts,
                 new_pidx_2 + num_single_cloth_verts);
  });

  const int base_vtx = num_single_cloth_verts * 2;
  const int base_faces = faces.rows() * 2;
  const scalar rot_angle = M_PI / (scalar)(num_cloth_edge_discretization);

  for (int iCloth = 0; iCloth < num_clothes; ++iCloth) {
    // compute boundary vector
    const std::vector<int> loop = m_face_loops[iCloth];
    const int num_bverts = loop.size();
    const int cloth_edge_vert_base =
        base_vtx + lv_base[iCloth] * (num_cloth_edge_discretization - 1);
    const int cloth_edge_face_base =
        base_faces + lv_base[iCloth] * num_cloth_edge_discretization * 2;

    threadutils::for_each(0, num_bverts, [&](int idx) {
      const int pidx = loop[idx];
      const int prev_pidx = loop[mathutils::mod_floor(idx - 1, num_bverts)];
      const int next_idx = mathutils::mod_floor(idx + 1, num_bverts);
      const int next_pidx = loop[mathutils::mod_floor(idx + 1, num_bverts)];

      const Vector3s e =
          (x.segment<3>(next_pidx * 4) - x.segment<3>(prev_pidx * 4))
              .normalized();

      for (int i = 1; i < num_cloth_edge_discretization; ++i) {
        Vector3s dir = Eigen::AngleAxis<scalar>(-rot_angle * (scalar)i, e) *
                       vert_norms[pidx];
        const scalar interp = (scalar)i / (scalar)num_cloth_edge_discretization;

        Vector3s pos = x.segment<3>(pidx * 4) +
                       dir * (radius(pidx * 2 + 0) * (1.0 - interp) +
                              radius(pidx * 2 + 1) * interp);
        const int new_pidx = cloth_edge_vert_base +
                             idx * (num_cloth_edge_discretization - 1) +
                             (i - 1);

        data->m_dbl_face_cloth_vertices[new_pidx] = pos;
        data->m_dbl_face_cloth_vertices_rest[new_pidx] =
            rest_x.segment<3>(pidx * 4);
        data->m_dbl_face_cloth_vertices_central[new_pidx] =
            x.segment<3>(pidx * 4);
        data->m_dbl_face_cloth_dir[new_pidx] = 1;
        data->m_dbl_face_cloth_group[new_pidx] = group[pidx];
        data->m_dbl_face_cloth_sat[new_pidx] =
            fluid_vol(pidx) / std::max(1e-16, vol(pidx));
      }

      for (int i = 0; i < num_cloth_edge_discretization; ++i) {
        int v_prev =
            (i == 0) ? (cloth_indicator[pidx] - 1)
                     : (cloth_edge_vert_base +
                        idx * (num_cloth_edge_discretization - 1) + (i - 1));
        int v_next = (i == num_cloth_edge_discretization - 1)
                         ? (num_single_cloth_verts + cloth_indicator[pidx] - 1)
                         : (cloth_edge_vert_base +
                            idx * (num_cloth_edge_discretization - 1) + i);

        int vn_prev =
            (i == 0)
                ? (cloth_indicator[next_pidx] - 1)
                : (cloth_edge_vert_base +
                   next_idx * (num_cloth_edge_discretization - 1) + (i - 1));
        int vn_next =
            (i == num_cloth_edge_discretization - 1)
                ? (num_single_cloth_verts + cloth_indicator[next_pidx] - 1)
                : (cloth_edge_vert_base +
                   next_idx * (num_cloth_edge_discretization - 1) + i);

        data->m_dbl_face_cloth_indices[cloth_edge_face_base +
                                       idx *
                                           (num_cloth_edge_discretization * 2) +
                                       i * 2 + 0] =
            Vector3i(v_prev, vn_next, v_next);
        data->m_dbl_face_cloth_indices[cloth_edge_face_base +
                                       idx *
                                           (num_cloth_edge_discretization * 2) +
                                       i * 2 + 1] =
            Vector3i(v_prev, vn_prev, vn_next);
      }
    });
  }
}
