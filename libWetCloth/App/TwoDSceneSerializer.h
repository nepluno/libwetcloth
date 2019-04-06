//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef __TWO_D_SCENE_SERIALIZER_H__
#define __TWO_D_SCENE_SERIALIZER_H__

#include <fstream>
#include <iostream>

#include "TwoDScene.h"
#include "StringUtilities.h"

struct SerializePosPacket
{
    std::string fn_pos;
    VectorXs m_pos;
    MatrixXs m_d_gauss;
};

struct SerializePacket
{
    std::string fn_clothes;
    std::string fn_hairs;
    std::string fn_fluid;
    std::string fn_internal_boundaries;
    std::string fn_external_boundaries;
    std::string fn_springs;
    
    std::vector< Vector3i > m_dbl_face_cloth_indices;
    std::vector< Vector3s > m_dbl_face_cloth_vertices;
    std::vector< Vector3s > m_dbl_face_cloth_vertices_rest;
    std::vector< Vector3s > m_dbl_face_cloth_vertices_central;
    std::vector< scalar > m_dbl_face_cloth_sat;
    std::vector< int > m_dbl_face_cloth_group;
    std::vector< int > m_dbl_face_cloth_dir;
    
    std::vector< Vector3s > m_hair_vertices;
    std::vector< Vector3s > m_hair_vertices_rest;
    std::vector< Vector2i > m_hair_indices;
    std::vector< Vector2s > m_hair_radii;
    std::vector< scalar > m_hair_sat;
    std::vector< int > m_hair_group;
    
    std::vector< Vector3s > m_attach_spring_vertices;
    
    std::vector< Vector3s > m_fluid_vertices;
    std::vector< scalar > m_fluid_radii;
    
    std::vector< Vector3i > m_external_indices;
    std::vector< Vector3s > m_external_vertices;
    
    std::vector< Vector3i > m_internal_indices;
    std::vector< Vector3s > m_internal_vertices;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class TwoDSceneSerializer
{
    std::vector< std::vector<int> > m_face_loops;
public:
  void serializeScene( TwoDScene& scene,
                      const std::string& fn_clothes,
                      const std::string& fn_hairs,
                      const std::string& fn_fluid,
                      const std::string& fn_internal_boundaries,
                      const std::string& fn_external_boundaries,
                      const std::string& fn_springs);
    
  void serializePositionOnly( TwoDScene& scene,
                             const std::string& fn_pos );

  void loadPosOnly( TwoDScene& scene, std::ifstream& inputstream );
    
    void initializeFaceLoops(const TwoDScene& scene);
    
    void updateDoubleFaceCloth(const TwoDScene& scene, SerializePacket* data);
    void updateHairs(const TwoDScene& scene, SerializePacket* data);
    void updateFluid(const TwoDScene& scene, SerializePacket* data);
    void updateMesh(const TwoDScene& scene, SerializePacket* data);
    void updateAttachSprings(const TwoDScene& scene, SerializePacket* data);
};

#endif
