//
// This file is part of the libWetCloth open source project
//
// The code is licensed solely for academic and non-commercial use under the
// terms of the Clear BSD License. The terms of the Clear BSD License are
// provided below. Other licenses may be obtained by contacting the faculty
// of the Columbia Computer Graphics Group or a Columbia University licensing officer.
//
// We would like to hear from you if you appreciate this work.
//
// The Clear BSD License
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//  list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
//  to endorse or promote products derived from this software without specific
//  prior written permission.
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

#include "CohesionForce.h"
#include "TwoDScene.h"

CohesionForce::CohesionForce( const std::shared_ptr<TwoDScene>& scene_ptr )
: Force()
, m_scene(scene_ptr)
{
}

CohesionForce::~CohesionForce()
{}

void CohesionForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
}

void CohesionForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
    const std::vector< std::vector<RayTriInfo> >& intersections = m_scene->getIntersections();
    
    const int num_edges = m_scene->getNumEdges();
    const int num_faces = m_scene->getNumFaces();
    const MatrixXi& faces = m_scene->getFaces();
    const MatrixXi& edges = m_scene->getEdges();
    const std::vector<int>& surfels = m_scene->getSurfels();
    const std::vector< Vector3s >& face_weights = m_scene->getFaceWeights();
    for(const std::vector<RayTriInfo>& gauss_info : intersections)
    {
        for(const RayTriInfo& info : gauss_info) {
            Vector3s w0, w1, x_start, x_end;
            Vector3i ele0, ele1;
            
            int num_elem_0, num_elem_1;
            
            if(info.start_geo_id < num_edges) {
                w0 = Vector3s(0.5, 0.5, 0.0);
                const Vector2iT& e = edges.row(info.start_geo_id);
                x_start = x.segment<3>(e(0) * 4) * w0(0) + x.segment<3>(e(1) * 4) * w0(1);
                ele0 = Vector3i(e(0), e(1), -1);
                num_elem_0 = 2;
            } else if(info.start_geo_id < num_edges + num_faces) {
                const int fidx_start = info.start_geo_id - num_edges;
                const Vector3iT& f = faces.row(fidx_start);
                w0 = face_weights[fidx_start];
                x_start = x.segment<3>(f(0) * 4) * w0(0) + x.segment<3>(f(1) * 4) * w0(1) + x.segment<3>(f(2) * 4) * w0(2);
                ele0 = f.transpose();
                num_elem_0 = 3;
            } else {
                continue;
            }
            
            if(info.intersect_geo_id < num_edges) {
                w1 = Vector3s(info.uv(0), 1.0 - info.uv(0), 0.0);
                const Vector2iT& e = edges.row(info.intersect_geo_id);
                x_end = x.segment<3>(e(0) * 4) * w1(0) + x.segment<3>(e(1) * 4) * w1(1);
                ele1 = Vector3i(e(0), e(1), -1);
                num_elem_1 = 2;
            } else if(info.intersect_geo_id < num_edges + num_faces) {
                const int fidx_end = info.intersect_geo_id - num_edges;
                const Vector3iT& f = faces.row(fidx_end);
                w1 = Vector3s(info.uv(0), info.uv(1), 1.0 - info.uv(0) - info.uv(1));
                x_end = x.segment<3>(f(0) * 4) * w1(0) + x.segment<3>(f(1) * 4) * w1(1) + x.segment<3>(f(2) * 4) * w1(2);
                ele1 = f.transpose();
                num_elem_1 = 3;
            } else {
                const int sidx_end = info.intersect_geo_id - num_edges - num_faces;
                w1 = Vector3s(1, 0, 0);
                x_end = x.segment<3>(surfels[sidx_end] * 4);
                ele1 = Vector3i(-1, -1, -1);
                num_elem_1 = 0;
            }
            
            Vector3s nhat = x_start - x_end;
            const scalar d = std::max(1e-12, nhat.norm());
            nhat /= d;
            
            const scalar dE = info.c1 * info.weight;
            const Vector3s gE = dE * nhat;
            
            for(int r = 0; r < num_elem_0; ++r) {
                gradE.segment<3>(ele0(r) * 4) += gE * w0(r);
            }
            
            for(int r = 0; r < num_elem_1; ++r) {
                gradE.segment<3>(ele1(r) * 4) += -gE * w1(r);
            }
        }
    }
}

void CohesionForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
    const std::vector< std::vector<RayTriInfo> >& intersections = m_scene->getIntersections();
    
    const int num_edges = m_scene->getNumEdges();
    const int num_faces = m_scene->getNumFaces();
    const MatrixXi& faces = m_scene->getFaces();
    const MatrixXi& edges = m_scene->getEdges();
    const std::vector<int>& surfels = m_scene->getSurfels();
    const std::vector< Vector3s >& face_weights = m_scene->getFaceWeights();
    const int num_gauss = intersections.size();
    
    threadutils::for_each(0, num_gauss, [&] (int gidx) {
        const std::vector<RayTriInfo>& gauss_info = intersections[gidx];
        const std::vector<int>& offsets = m_hess_offsets[gidx];
        
        int k = 0;
        for(const RayTriInfo& info : gauss_info) {
            Vector3s w0, w1, x_start, x_end;
            Vector3i ele0, ele1;
            
            int num_elem_0, num_elem_1;
            
            if(info.start_geo_id < num_edges) {
                w0 = Vector3s(0.5, 0.5, 0.0);
                const Vector2iT& e = edges.row(info.start_geo_id);
                x_start = x.segment<3>(e(0) * 4) * w0(0) + x.segment<3>(e(1) * 4) * w0(1);
                ele0 = Vector3i(e(0), e(1), -1);
                num_elem_0 = 2;
            } else if(info.start_geo_id < num_edges + num_faces) {
                const int fidx_start = info.start_geo_id - num_edges;
                const Vector3iT& f = faces.row(fidx_start);
                w0 = face_weights[fidx_start];
                x_start = x.segment<3>(f(0) * 4) * w0(0) + x.segment<3>(f(1) * 4) * w0(1) + x.segment<3>(f(2) * 4) * w0(2);
                ele0 = f.transpose();
                num_elem_0 = 3;
            } else {
                continue;
            }
            
            if(info.intersect_geo_id < num_edges) {
                w1 = Vector3s(info.uv(0), 1.0 - info.uv(0), 0.0);
                const Vector2iT& e = edges.row(info.intersect_geo_id);
                x_end = x.segment<3>(e(0) * 4) * w1(0) + x.segment<3>(e(1) * 4) * w1(1);
                ele1 = Vector3i(e(0), e(1), -1);
                num_elem_1 = 2;
            } else if(info.intersect_geo_id < num_edges + num_faces) {
                const int fidx_end = info.intersect_geo_id - num_edges;
                const Vector3iT& f = faces.row(fidx_end);
                w1 = Vector3s(info.uv(0), info.uv(1), 1.0 - info.uv(0) - info.uv(1));
                x_end = x.segment<3>(f(0) * 4) * w1(0) + x.segment<3>(f(1) * 4) * w1(1) + x.segment<3>(f(2) * 4) * w1(2);
                ele1 = f.transpose();
                num_elem_1 = 3;
            } else {
                const int sidx_end = info.intersect_geo_id - num_edges - num_faces;
                w1 = Vector3s(1, 0, 0);
                x_end = x.segment<3>(surfels[sidx_end] * 4);
                ele1 = Vector3i(-1, -1, -1);
                num_elem_1 = 0;
            }
            
            Vector3s nhat = x_start - x_end;
            const scalar d = std::max(1e-12, nhat.norm());
            
            nhat /= d;
            
            Matrix3s J = (info.c1 * (Matrix3s::Identity() - nhat * nhat.transpose()) / d) * info.weight;
            
            int base_idx = offsets[k];
            
            for(int r = 0; r < num_elem_0; ++r) for(int s = 0; s < num_elem_0; ++s)
            {
                int element_idx = r * num_elem_0 + s;
                for(int i = 0; i < 3; ++i) for(int j = 0; j < 3; ++j)
                {
                    hessE[hessE_index + base_idx + element_idx * 9 + i * 3 + j] = Triplets(4 * ele0(r) + i, 4 * ele0(s) + j, J(i, j) * w0(r) * w0(s));
                }
            }
            
            base_idx += num_elem_0 * num_elem_0 * 9;
            
            for(int r = 0; r < num_elem_0; ++r) for(int s = 0; s < num_elem_1; ++s)
            {
                int element_idx = r * num_elem_1 + s;
                for(int i = 0; i < 3; ++i) for(int j = 0; j < 3; ++j)
                {
                    hessE[hessE_index + base_idx + element_idx * 9 + i * 3 + j] = Triplets(4 * ele0(r) + i, 4 * ele1(s) + j, J(i, j) * w0(r) * w1(s));
                }
            }
            
            base_idx += num_elem_0 * num_elem_1 * 9;
            
            for(int r = 0; r < num_elem_1; ++r) for(int s = 0; s < num_elem_0; ++s)
            {
                int element_idx = r * num_elem_0 + s;
                for(int i = 0; i < 3; ++i) for(int j = 0; j < 3; ++j)
                {
                    hessE[hessE_index + base_idx + element_idx * 9 + i * 3 + j] = Triplets(4 * ele1(r) + i, 4 * ele0(s) + j, J(i, j) * w1(r) * w0(s));
                }
            }
            
            base_idx += num_elem_1 * num_elem_0 * 9;
            
            for(int r = 0; r < num_elem_1; ++r) for(int s = 0; s < num_elem_1; ++s)
            {
                int element_idx = r * num_elem_1 + s;
                for(int i = 0; i < 3; ++i) for(int j = 0; j < 3; ++j)
                {
                    hessE[hessE_index + base_idx + element_idx * 9 + i * 3 + j] = Triplets(4 * ele1(r) + i, 4 * ele1(s) + j, J(i, j) * w1(r) * w1(s));
                }
            }
            ++k;
        }
    });
}

void CohesionForce::preCompute()
{
}

void CohesionForce::updateStartState()
{
}

int CohesionForce::numHessX()
{
    const int num_edges = m_scene->getNumEdges();
    const int num_faces = m_scene->getNumFaces();
    
    const std::vector< std::vector<RayTriInfo> >& intersections = m_scene->getIntersections();
    m_hess_offsets.resize(intersections.size());
    
    int num_inters = 0;
    int k = 0;
    for(const std::vector<RayTriInfo>& gauss_info : intersections)
    {
        const int num_info = gauss_info.size();
        m_hess_offsets[k].resize(0);
        for(int i = 0; i < num_info; ++i)
        {
            int s0, s1;
            if(gauss_info[i].start_geo_id < num_edges) s0 = 2;
            else if(gauss_info[i].start_geo_id < num_edges + num_faces) s0 = 3;
            else s0 = 0;
            
            if(gauss_info[i].intersect_geo_id < num_edges) s1 = 2;
            else if(gauss_info[i].intersect_geo_id < num_edges + num_faces) s1 = 3;
            else s1 = 0;
            
            const int num_hess = (s0 + s1) * (s0 + s1) * 9;
            
            m_hess_offsets[k].push_back(num_inters);
            num_inters += num_hess;
        }
        k++;
    }
    
	return num_inters;
}

Force* CohesionForce::createNewCopy()
{
	return new CohesionForce(*this);
}

int CohesionForce::flag() const
{
	return 1;
}
