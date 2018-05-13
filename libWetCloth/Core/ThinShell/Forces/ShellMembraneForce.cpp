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


#include "ShellMembraneForce.h"
#include "../../ThreadUtils.h"
#include <iostream>

ShellMembraneForce::~ShellMembraneForce()
{}

void ShellMembraneForce::perFaceNormals( const VectorXs& V, const MatrixXi& F, MatrixXs& N )
{
	N.resize(F.rows(), 3);
	
	int Frows = F.rows();
	
	threadutils::for_each(0, Frows, [&] (int i) {
		const Vector3s& v1 = V.segment<3>(F(i, 1) * 4) - V.segment<3>(F(i, 0) * 4);
		const Vector3s& v2 = V.segment<3>(F(i, 2) * 4) - V.segment<3>(F(i, 0) * 4);
		
		N.row(i) = (v1.cross(v2)).transpose();
		const scalar r = N.row(i).norm();
		
		if(r > 0.0) {
			N.row(i) /= r;
		}
	});
}

ShellMembraneForce::ShellMembraneForce(const VectorXs & rest_pos,
									   const VectorXs & pos,
									   const MatrixXi & F,
									   const VectorXs & triangle_rest_area,
									   const scalar& young_modulus,
									   const scalar& viscous_modulus,
									   const scalar& poisson_ratio)
: m_rest_pos(rest_pos),
m_pos(pos),
m_F(F),
m_triangle_rest_area(triangle_rest_area),
m_young_modulus(young_modulus),
m_viscous_modulus(viscous_modulus),
m_poisson_ratio(poisson_ratio)
{
	perFaceNormals(m_rest_pos, m_F, m_triangle_normals);
	m_membrane_ru.resize(m_F.rows(), 3);
	m_membrane_rv.resize(m_F.rows(), 3);
	
	threadutils::for_each(0, (int) m_F.rows(), [&] (int f) {
		const Vector3s& x0 = m_rest_pos.segment<3>(m_F(f,0) * 4).transpose();
		const Vector3s& x1 = m_rest_pos.segment<3>(m_F(f,1) * 4).transpose();
		const Vector3s& x2 = m_rest_pos.segment<3>(m_F(f,2) * 4).transpose();
		
		const Vector3s normal = m_triangle_normals.row(f);
		
		//Define (u,v) coords by LOCALLY rotating the triangle to align it with the Z axis
		Eigen::Quaternion<scalar> rot = Eigen::Quaternion<scalar>::FromTwoVectors(normal, Vector3s::UnitZ());
		
		// compute uv of each vertex
		Vector2s uvi = Vector2s::Zero();
		Vector2s uvj = (rot * (x1 - x0)).segment(0,2);
		Vector2s uvk = (rot * (x2 - x0)).segment(0,2);
		
		//Determine vertex weights for strain computation
		const scalar dinv = 1./(uvi(0) * (uvj(1) - uvk(1)) + uvj(0) * (uvk(1) - uvi(1)) + uvk(0) * (uvi(1) - uvj(1)));
		m_membrane_ru.row(f) << dinv * (uvj(1) - uvk(1)), dinv * (uvk(1) - uvi(1)), dinv * (uvi(1) - uvj(1));
		m_membrane_rv.row(f) << dinv * (uvk(0) - uvj(0)), dinv * (uvi(0) - uvk(0)), dinv * (uvj(0) - uvi(0));
	});
	
	m_membrane_material_tensor << 1, m_poisson_ratio, 0, m_poisson_ratio, 1, 0, 0, 0, 0.5 * (1 - m_poisson_ratio);
	m_membrane_material_viscous_tensor = m_membrane_material_tensor * (m_viscous_modulus / (1 - m_poisson_ratio * m_poisson_ratio));
	m_membrane_material_tensor = m_membrane_material_tensor * (m_young_modulus / (1 - m_poisson_ratio * m_poisson_ratio));
}

void ShellMembraneForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	std::cerr << "NOT IMPLEMENTED!" << std::endl;
}

void ShellMembraneForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
	for (int f = 0; f < m_F.rows(); ++f) {
		const scalar psi_base = (psi(m_F(f,0)) + psi(m_F(f,1)) + psi(m_F(f,2))) / 3.0;
		const scalar psi_coeff = pow(psi_base, lambda);
		
		const Vector3s& x0 = x.segment<3>(m_F(f,0) * 4);
		const Vector3s& x1 = x.segment<3>(m_F(f,1) * 4);
		const Vector3s& x2 = x.segment<3>(m_F(f,2) * 4);
		
		const Vector3s U = x0*m_membrane_ru(f,0) + x1*m_membrane_ru(f,1) + x2*m_membrane_ru(f,2);
		const Vector3s V = x0*m_membrane_rv(f,0) + x1*m_membrane_rv(f,1) + x2*m_membrane_rv(f,2);
		
		const Vector3s strain = Vector3s(0.5*(U.dot(U) - 1), 0.5*(V.dot(V) - 1), U.dot(V));
		const Vector3s stress = m_membrane_material_tensor*strain;
		
		for(int i = 0; i < 3; ++i) {
			gradE.segment<3>(m_F(f, i) * 4) += psi_coeff *
			m_triangle_rest_area(f) * (stress(0) * (m_membrane_ru(f,i) * U) +
									   stress(1) * (m_membrane_rv(f,i) * V) +
									   stress(2) * (m_membrane_ru(f,i) * V +
													m_membrane_rv(f,i) * U));
		}
		
		if(m_viscous_modulus > 0.0) {
			const Vector3s& sx0 = m_start_pos.segment<3>(m_F(f,0) * 4);
			const Vector3s& sx1 = m_start_pos.segment<3>(m_F(f,1) * 4);
			const Vector3s& sx2 = m_start_pos.segment<3>(m_F(f,2) * 4);
			
			const Vector3s sU = sx0*m_membrane_ru(f,0) + sx1*m_membrane_ru(f,1) + sx2*m_membrane_ru(f,2);
			const Vector3s sV = sx0*m_membrane_rv(f,0) + sx1*m_membrane_rv(f,1) + sx2*m_membrane_rv(f,2);
			
			const Vector3s viscous_strain = Vector3s(0.5*(U.dot(U) - sU.dot(sU)), 0.5*(V.dot(V) - sV.dot(sV)), U.dot(V) - sU.dot(sV));
			const Vector3s viscous_stress = m_membrane_material_viscous_tensor*strain;
			
			for(int i = 0; i < 3; ++i) {
				gradE.segment<3>(m_F(f, i) * 4) += psi_coeff *
				m_triangle_rest_area(f) * (viscous_stress(0) * (m_membrane_ru(f,i) * U) +
										   viscous_stress(1) * (m_membrane_rv(f,i) * V) +
										   viscous_stress(2) * (m_membrane_ru(f,i) * V +
														m_membrane_rv(f,i) * U));
			}
		}
	}
}

void ShellMembraneForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	threadutils::for_each(0, (int) m_F.rows(), [&] (int f) {
		const scalar psi_base = (psi(m_F(f,0)) + psi(m_F(f,1)) + psi(m_F(f,2))) / 3.0;
		const scalar psi_coeff = pow(psi_base, lambda);
		const Vector3s& x0 = x.segment<3>(m_F(f, 0) * 4);
		const Vector3s& x1 = x.segment<3>(m_F(f, 1) * 4);
		const Vector3s& x2 = x.segment<3>(m_F(f, 2) * 4);
		
		Vector3s U = x0*m_membrane_ru(f,0) + x1*m_membrane_ru(f,1) + x2*m_membrane_ru(f,2);
		Vector3s V = x0*m_membrane_rv(f,0) + x1*m_membrane_rv(f,1) + x2*m_membrane_rv(f,2);
		
		Vector3s strain = Vector3s(0.5 * (U.dot(U) - 1), 0.5 * (V.dot(V) - 1), U.dot(V));
		Vector3s stress = m_membrane_material_tensor*strain;
		
		const Vector3s& sx0 = m_start_pos.segment<3>(m_F(f,0) * 4);
		const Vector3s& sx1 = m_start_pos.segment<3>(m_F(f,1) * 4);
		const Vector3s& sx2 = m_start_pos.segment<3>(m_F(f,2) * 4);
		
		const Vector3s sU = sx0*m_membrane_ru(f,0) + sx1*m_membrane_ru(f,1) + sx2*m_membrane_ru(f,2);
		const Vector3s sV = sx0*m_membrane_rv(f,0) + sx1*m_membrane_rv(f,1) + sx2*m_membrane_rv(f,2);
		
		const Vector3s viscous_strain = Vector3s(0.5*(U.dot(U) - sU.dot(sU)), 0.5*(V.dot(V) - sV.dot(sV)), U.dot(V) - sU.dot(sV));
		const Vector3s viscous_stress = m_membrane_material_viscous_tensor*strain;
		
		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				Matrix3s dphidx_membrane = (psi_coeff * m_triangle_rest_area(f))
				*(m_membrane_material_tensor(0,0)*(m_membrane_ru(f,i)*m_membrane_ru(f,j)*U*U.transpose())
				+ m_membrane_material_tensor(1,1)*(m_membrane_rv(f,i)*m_membrane_rv(f,j)*V*V.transpose())
				+ m_membrane_material_tensor(2,2)*( (m_membrane_rv(f,i)*m_membrane_ru(f,j)*U*V.transpose()
				+ m_membrane_ru(f,i)*m_membrane_rv(f,j)*V*U.transpose())
				+ (m_membrane_rv(f,i)*m_membrane_rv(f,j)*U*U.transpose() + m_membrane_ru(f,i)*m_membrane_ru(f,j)*V*V.transpose()) )
				+ m_membrane_material_tensor(0,1)*(m_membrane_ru(f,i)*m_membrane_rv(f,j)*U*V.transpose())
				+ m_membrane_material_tensor(1,0)*(m_membrane_rv(f,i)*m_membrane_ru(f,j)*V*U.transpose())
				+ (stress(0)*(m_membrane_ru(f,i)*m_membrane_ru(f,j)) + stress(1)*(m_membrane_rv(f,i)*m_membrane_rv(f,j)) + stress(2)*(m_membrane_ru(f,i)*m_membrane_rv(f,j) + m_membrane_rv(f,i)*m_membrane_ru(f,j)))*Matrix3s::Identity());
				
				if(m_viscous_modulus > 0.0) {
					dphidx_membrane += (psi_coeff * m_triangle_rest_area(f))
					  *(m_membrane_material_viscous_tensor(0,0)*(m_membrane_ru(f,i)*m_membrane_ru(f,j)*U*U.transpose())
					  + m_membrane_material_viscous_tensor(1,1)*(m_membrane_rv(f,i)*m_membrane_rv(f,j)*V*V.transpose())
					  + m_membrane_material_viscous_tensor(2,2)*((m_membrane_rv(f,i)*m_membrane_ru(f,j)*U*V.transpose()
														   + m_membrane_ru(f,i)*m_membrane_rv(f,j)*V*U.transpose())
														 + (m_membrane_rv(f,i)*m_membrane_rv(f,j)*U*U.transpose() + m_membrane_ru(f,i)*m_membrane_ru(f,j)*V*V.transpose()) )
					  + m_membrane_material_viscous_tensor(0,1)*(m_membrane_ru(f,i)*m_membrane_rv(f,j)*U*V.transpose())
					  + m_membrane_material_viscous_tensor(1,0)*(m_membrane_rv(f,i)*m_membrane_ru(f,j)*V*U.transpose())
					  + (viscous_stress(0)*(m_membrane_ru(f,i)*m_membrane_ru(f,j)) + viscous_stress(1)*(m_membrane_rv(f,i)*m_membrane_rv(f,j)) + viscous_stress(2)*(m_membrane_ru(f,i)*m_membrane_rv(f,j) + m_membrane_rv(f,i)*m_membrane_ru(f,j)))*Matrix3s::Identity());
				}
				
				int base_idx = f * 81 + (i * 3 + j) * 9;
				
				for(int s = 0; s < 3; ++s) {
					for(int r = 0; r < 3; ++r) {
						hessE[hessE_index + base_idx + s * 3 + r] = Triplets(m_F(f, i) * 4 + r, m_F(f, j) * 4 + s, dphidx_membrane(r, s));
					}
				}
			}
		}
	});
}

int ShellMembraneForce::numHessX()
{
	return m_F.rows() * 81;
}

void ShellMembraneForce::preCompute()
{
	perFaceNormals(m_pos, m_F, m_triangle_normals);
}

void ShellMembraneForce::updateStartState()
{
	m_start_pos = m_pos;
}

Force* ShellMembraneForce::createNewCopy()
{
	return new ShellMembraneForce(*this);
}

int ShellMembraneForce::flag() const
{
	return 1;
}
