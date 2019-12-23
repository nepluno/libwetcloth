//
// This file is part of the libWetCloth open source project
//
// Copyright 2016 Gabriel Cirio
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#include "ThinShellForce.h"
#include "../DER/StrandParameters.h"

#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unique_edge_map.h>
#include <igl/edge_flaps.h>
#include <igl/per_face_normals.h>
#include <igl/qslim.h>
#include <igl/decimate.h>

#include <memory>

bool ThinShellForce::parallelized() const
{
	return true;
}

ThinShellForce::~ThinShellForce()
{}

ThinShellForce::ThinShellForce(const std::shared_ptr<TwoDScene>& scene, const std::vector< Vector3i >& faces, const int& parameterIndex, int globalIndex)
: m_scene(scene)
{
	const VectorXs& rest_pos = scene->getRestPos();
	const int num_particles = scene->getNumParticles();
	
	m_F.resize(faces.size(), 3);
	m_triangle_rest_areas.resize(faces.size());
	const int num_faces = (int) faces.size();
	for(int i = 0; i < num_faces; ++i)  {
		m_F.row(i) = faces[i].transpose();
		
		const Vector3s x0 = rest_pos.segment<3>(m_F(i, 0) * 4);
		const Vector3s x1 = rest_pos.segment<3>(m_F(i, 1) * 4);
		const Vector3s x2 = rest_pos.segment<3>(m_F(i, 2) * 4);
		
		double area = 0.5*(x1 - x0).cross(x2 - x0).norm();
		m_triangle_rest_areas(i) = area;
	}
	
	igl::vertex_triangle_adjacency(num_particles, m_F, m_per_node_triangles, m_per_node_triangles_node_local_index);
	igl::triangle_triangle_adjacency(m_F, m_per_triangle_triangles, m_per_triangle_triangles_edge_local_index);
	igl::unique_edge_map(m_F, m_E_directed, m_E_unique, m_map_edge_directed_to_unique, m_map_edge_unique_to_directed);
	igl::edge_flaps(m_F, m_E_unique, m_map_edge_directed_to_unique, m_per_unique_edge_triangles, m_per_unique_edge_triangles_local_corners);
	
	m_per_node_edges.clear();
	m_node_neighbors.clear();
	
	m_per_node_edges.resize(num_particles);
	m_node_neighbors.resize(num_particles);
	
	for (int i = 0; i < m_E_unique.rows(); i++)
	{
		m_per_node_edges[m_E_unique(i,0)].push_back(i);
		m_per_node_edges[m_E_unique(i,1)].push_back(i);
	}
	
	for (int i = 0; i < num_particles; i++)
	{
		const int num_ne = (int) m_per_node_edges[i].size();
		for (int j=0; j < num_ne; j++)
		{
			int node0 = m_E_unique(m_per_node_edges[i][j], 0);
			int node1 = m_E_unique(m_per_node_edges[i][j], 1);
			
			if (node0 == i) m_node_neighbors[node0].push_back(node1);
			else m_node_neighbors[node1].push_back(node0);
		}
	}
	
	m_per_triangles_unique_edges.resize(m_F.rows(), 3);
	
	for (unsigned int e=0; e<m_per_unique_edge_triangles.rows(); e++)
	{
		int va = m_E_unique(e,0);
		int vb = m_E_unique(e,1);
		
		int f = m_per_unique_edge_triangles(e,0);
		if (f != -1)
		{
			if ( ( (m_F(f,1) == va) && (m_F(f,2) == vb) ) || ( (m_F(f,2) == va) && (m_F(f,1) == vb) ) )
				m_per_triangles_unique_edges(f,0) = e;
			else if ( ( (m_F(f,0) == va) && (m_F(f,2) == vb) ) || ( (m_F(f,2) == va) && (m_F(f,0) == vb) ) )
				m_per_triangles_unique_edges(f,1) = e;
			else if ( ( (m_F(f,0) == va) && (m_F(f,1) == vb) ) || ( (m_F(f,1) == va) && (m_F(f,0) == vb) ) )
				m_per_triangles_unique_edges(f,2) = e;
			else
				std::cout << "IMPOSSIBLE!!!" << std::endl;
		}
		
		f = m_per_unique_edge_triangles(e,1);
		if (f != -1)
		{
			if ( ( (m_F(f,1) == va) && (m_F(f,2) == vb) ) || ( (m_F(f,2) == va) && (m_F(f,1) == vb) ) )
				m_per_triangles_unique_edges(f,0) = e;
			else if ( ( (m_F(f,0) == va) && (m_F(f,2) == vb) ) || ( (m_F(f,2) == va) && (m_F(f,0) == vb) ) )
				m_per_triangles_unique_edges(f,1) = e;
			else if ( ( (m_F(f,0) == va) && (m_F(f,1) == vb) ) || ( (m_F(f,1) == va) && (m_F(f,0) == vb) ) )
				m_per_triangles_unique_edges(f,2) = e;
			else
				std::cout << "IMPOSSIBLE!!!" << std::endl;
		}
	}
	
	auto& params = scene->getStrandParameters(parameterIndex);
	
	const scalar poisson_ratio = params->m_youngsModulus.get() / (2.0 * params->m_shearModulus.get()) - 1.0;
	
	m_forces.push_back( std::make_shared<ShellMembraneForce>( 
		rest_pos, scene->getX(), 
		m_F, m_triangle_rest_areas, 
		params->m_youngsModulus.get(),
		params->m_viscousBendingCoefficientBase,
		poisson_ratio,
		params->m_physicalRadius.get()(0) + params->m_physicalRadius.get()(1),
		params->m_accumulateWithViscous && !params->m_accumulateViscousOnlyForBendingModes) );
	m_forces.push_back( std::make_shared<ShellBendingForce>(
		scene->getX(), rest_pos, 
		m_F, m_triangle_rest_areas, 
		m_E_unique, m_per_unique_edge_triangles, 
		m_per_unique_edge_triangles_local_corners, 
		m_per_triangles_unique_edges,
		params->m_youngsModulus.get(),
		params->m_viscousBendingCoefficientBase,
		poisson_ratio, 
		params->m_physicalRadius.get()(0) + params->m_physicalRadius.get()(1), 
		m_scene->getLiquidInfo().bending_scheme,
		params->m_accumulateWithViscous) );
}

void ThinShellForce::updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt )
{
	for(auto& force : m_forces)
	{
		force->updateMultipliers(x, vplus, m, psi, lambda, dt);
	}	
}


void ThinShellForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	for(auto& force : m_forces)
	{
		force->addEnergyToTotal(x, v, m, psi, lambda, E);
	}
}

void ThinShellForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
	for(auto& force : m_forces)
	{
		force->addGradEToTotal(x, v, m, psi, lambda, gradE);
	}
}

void ThinShellForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	int idx = hessE_index;
	for(auto& force : m_forces)
	{
		force->addHessXToTotal(x, v, m, psi, lambda, hessE, idx, dt);
		idx += force->numHessX();
	}
}

int ThinShellForce::numHessX()
{
	int idx = 0;
	for(auto& force : m_forces)
	{
		idx += force->numHessX();
	}
	return idx;
}

void ThinShellForce::preCompute()
{
	for(auto& force : m_forces)
	{
		force->preCompute();
	}
}

void ThinShellForce::updateStartState()
{
	for(auto& force : m_forces)
	{
		force->updateStartState();
	}
}

Force* ThinShellForce::createNewCopy()
{
	return new ThinShellForce(*this);
}

int ThinShellForce::flag() const
{
	return 1;
}

