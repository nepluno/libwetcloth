//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TWO_D_SCENE_SERIALIZER_H
#define TWO_D_SCENE_SERIALIZER_H

#include <fstream>
#include <iostream>

#include "StringUtilities.h"
#include "TwoDScene.h"

struct SerializePosPacket {
  std::string fn_pos;
  VectorXs m_pos;
  MatrixXs m_d_gauss;
};

struct SerializePacket {
  std::string fn_clothes;
  std::string fn_hairs;
  std::string fn_fluid;
  std::string fn_internal_boundaries;
  std::string fn_external_boundaries;
  std::string fn_springs;

  std::vector<Vector3i> m_dbl_face_cloth_indices;
  std::vector<Vector3s> m_dbl_face_cloth_vertices;
  std::vector<Vector3s> m_dbl_face_cloth_vertices_rest;
  std::vector<Vector3s> m_dbl_face_cloth_vertices_central;
  std::vector<scalar> m_dbl_face_cloth_sat;
  std::vector<int> m_dbl_face_cloth_group;
  std::vector<int> m_dbl_face_cloth_dir;

  std::vector<Vector3s> m_hair_vertices;
  std::vector<Vector3s> m_hair_vertices_rest;
  std::vector<Vector2i> m_hair_indices;
  std::vector<Vector2s> m_hair_radii;
  std::vector<scalar> m_hair_sat;
  std::vector<int> m_hair_group;

  std::vector<Vector3s> m_attach_spring_vertices;

  std::vector<Vector3s> m_fluid_vertices;
  std::vector<scalar> m_fluid_radii;

  std::vector<Vector3i> m_external_indices;
  std::vector<Vector3s> m_external_vertices;

  std::vector<Vector3i> m_internal_indices;
  std::vector<Vector3s> m_internal_vertices;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class TwoDSceneSerializer {
  std::vector<std::vector<int> > m_face_loops;

 public:
  void serializeScene(TwoDScene& scene, const std::string& fn_clothes,
                      const std::string& fn_hairs, const std::string& fn_fluid,
                      const std::string& fn_internal_boundaries,
                      const std::string& fn_external_boundaries,
                      const std::string& fn_springs);

  void serializePositionOnly(TwoDScene& scene, const std::string& fn_pos);

  void loadPosOnly(TwoDScene& scene, std::ifstream& inputstream);

  void initializeFaceLoops(const TwoDScene& scene);

  void updateDoubleFaceCloth(const TwoDScene& scene, SerializePacket* data);
  void updateHairs(const TwoDScene& scene, SerializePacket* data);
  void updateFluid(const TwoDScene& scene, SerializePacket* data);
  void updateMesh(const TwoDScene& scene, SerializePacket* data);
  void updateAttachSprings(const TwoDScene& scene, SerializePacket* data);
};

#endif
