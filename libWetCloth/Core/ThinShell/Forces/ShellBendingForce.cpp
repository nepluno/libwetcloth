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


#include "ShellBendingForce.h"
#include <iostream>
#include "../../ThreadUtils.h"
#undef isnan
#undef isinf

ShellBendingForce::~ShellBendingForce()
{}

ShellBendingForce::ShellBendingForce(const VectorXs & pos,
                                     const VectorXs & rest_pos,
                                     const MatrixXi & F,
                                     const VectorXs & triangle_rest_area,
                                     const MatrixXi & E_unique,
                                     const MatrixXi & per_unique_edge_triangles,
                                     const MatrixXi & per_unique_edge_triangles_local_corners,
                                     const MatrixXi & per_triangles_unique_edges,
                                     const scalar& young_modulus,
                                     const scalar& viscous_modulus,
                                     const scalar& poisson_ratio,
                                     const scalar& thickness,
                                     int bending_mode)
: m_pos(pos), m_rest_pos(rest_pos), m_F(F), m_triangle_rest_area(triangle_rest_area), m_E_unique(E_unique),
m_per_unique_edge_triangles(per_unique_edge_triangles), m_per_unique_edge_triangles_local_corners(per_unique_edge_triangles_local_corners), m_per_triangles_unique_edges(per_triangles_unique_edges), m_young_modulus(young_modulus), m_viscous_modulus(viscous_modulus),
m_poisson_ratio(poisson_ratio), m_thickness(thickness), m_bending_mode(bending_mode)
{
	m_per_edge_rest_phi.resize(m_E_unique.rows());
	m_per_edge_rest_phi.setZero();

	m_bending_stiffness = m_young_modulus * m_thickness * m_thickness * m_thickness / (24. * (1. - m_poisson_ratio * m_poisson_ratio));
	m_viscous_stiffness = m_viscous_modulus * m_thickness * m_thickness * m_thickness / (24. * (1. - m_poisson_ratio * m_poisson_ratio));

	for(int e = 0; e < E_unique.rows(); ++e)
	{
		if ( (m_per_unique_edge_triangles(e, 0) == -1) || (m_per_unique_edge_triangles(e, 1) == -1) )
			continue;
		
		m_unique_edge_usable.push_back(e);
	}
	
	computeBendingRestPhi(m_rest_pos, m_per_edge_rest_phi);
}

void ShellBendingForce::computeBendingRestPhi(const VectorXs & rest_pos, VectorXs & rest_phis)
{
	auto l_bending_rest_phi = [&] (int idx[4], scalar & rest_phi)
	{
		auto l_compute_tantheta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, scalar & rest_phi)
		{
			scalar tan_half_theta = (n - n_tilde).norm()/(n + n_tilde).norm();
			scalar sign_angle = n.cross(n_tilde).dot(e0);
			if (std::isnan(sign_angle)) sign_angle = 1.;
			else sign_angle = sign_angle > 0 ? 1. : -1.;
			rest_phi = 2.*sign_angle*tan_half_theta;
		};
		
		auto l_compute_sintheta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, scalar & rest_phi)
		{
			scalar theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));
			rest_phi = std::sin(theta/2.);
		};
		
		auto l_compute_theta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, scalar & rest_phi)
		{
			scalar theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));
			rest_phi = theta;
		};
		
		const Vector3s& x0 = rest_pos.segment<3>(idx[0] * 4);
		const Vector3s& x1 = rest_pos.segment<3>(idx[1] * 4);
		const Vector3s& x2 = rest_pos.segment<3>(idx[2] * 4);
		const Vector3s& x3 = rest_pos.segment<3>(idx[3] * 4);
		
		Vector3s e0 = x2 - x1;
		Vector3s e1 = x0 - x2;
		Vector3s e1_tilde = x3 - x2;
		
		Vector3s n = e0.cross(e1);
		scalar n_length = n.norm();
		Vector3s n_tilde = -e0.cross(e1_tilde);
		scalar n_tilde_length = n_tilde.norm();
		n /= n_length;
		n_tilde /= n_tilde_length;
		
        switch (m_bending_mode) {
            case 0:
                l_compute_tantheta(n, n_tilde, e0, rest_phi);
                break;
            case 1:
                l_compute_sintheta(n, n_tilde, e0, rest_phi);
                break;
            default:
                l_compute_theta(n, n_tilde, e0, rest_phi);
                break;
        }
	};
	
	threadutils::for_each(0, (int) m_unique_edge_usable.size(), [&] (int i) {
		const int e = m_unique_edge_usable[i];
		
		assert ( (m_per_unique_edge_triangles(e,0) != -1) && (m_per_unique_edge_triangles(e,1) != -1) );
		
		int idx[4];
		idx[0] = m_F(m_per_unique_edge_triangles(e,0), m_per_unique_edge_triangles_local_corners(e,0));
		idx[1] = m_E_unique(e,0);
		idx[2] = m_E_unique(e,1);
		idx[3] = m_F(m_per_unique_edge_triangles(e,1), m_per_unique_edge_triangles_local_corners(e,1));
		
		l_bending_rest_phi(idx, rest_phis(e));
	});
}

void ShellBendingForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	std::cerr << "NOT IMPLEMENTED!" << std::endl;
}

void ShellBendingForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
	for (int e : m_unique_edge_usable)
	{
		assert ( (m_per_unique_edge_triangles(e,0) != -1) && (m_per_unique_edge_triangles(e,1) != -1) );
		
		int idx[4];
		idx[0] = m_F(m_per_unique_edge_triangles(e,0), m_per_unique_edge_triangles_local_corners(e,0));
		idx[1] = m_E_unique(e,0);
		idx[2] = m_E_unique(e,1);
		idx[3] = m_F(m_per_unique_edge_triangles(e,1), m_per_unique_edge_triangles_local_corners(e,1));
		
		const scalar psi_coeff = pow((psi(idx[0]) + psi(idx[1]) + psi(idx[2]) + psi(idx[3])) * 0.25, lambda);
		
		auto l_compute_dPsi_tantheta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, const scalar rest_phi, const scalar ka, const scalar start_phi, const scalar kb, scalar & dPsi_dTheta)
		{
			scalar tan_half_theta = (n - n_tilde).norm()/(n + n_tilde).norm();
			scalar sign_angle = n.cross(n_tilde).dot(e0);
			if (std::isnan(sign_angle)) sign_angle = 1.;
			else sign_angle = sign_angle > 0 ? 1. : -1.;
			scalar phi_theta = 2.*sign_angle*tan_half_theta;
			dPsi_dTheta = 2.*(ka*(phi_theta - rest_phi) + kb*(phi_theta - start_phi))*(1.+tan_half_theta*tan_half_theta);
		};
		
		auto l_compute_dPsi_sintheta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, const scalar rest_phi, const scalar ka, const scalar start_phi, const scalar kb, scalar & dPsi_dTheta)
		{
			scalar theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));
			scalar phi_theta = std::sin(theta/2.);
			dPsi_dTheta = 2.*(ka*(phi_theta - rest_phi) + kb*(phi_theta - start_phi))*(0.5*std::cos(theta/2.));
		};
		
		auto l_compute_dPsi_theta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, const scalar rest_phi, const scalar ka, const scalar start_phi, const scalar kb, scalar & dPsi_dTheta)
		{
			scalar theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));
			scalar phi_theta = theta;
			dPsi_dTheta = 2.*(ka*(phi_theta - rest_phi) + kb*(phi_theta - start_phi));
		};
		
		const Vector3s& x0 = x.segment<3>(idx[0] * 4);
		const Vector3s& x1 = x.segment<3>(idx[1] * 4);
		const Vector3s& x2 = x.segment<3>(idx[2] * 4);
		const Vector3s& x3 = x.segment<3>(idx[3] * 4);
		
		scalar restareas = m_triangle_rest_area(m_per_unique_edge_triangles(e,0)) + m_triangle_rest_area(m_per_unique_edge_triangles(e,1));
		scalar e0_rest_sqnorm = (m_rest_pos.segment<3>(idx[2] * 4) - m_rest_pos.segment<3>(idx[1] * 4)).squaredNorm();
		
		Vector3s e0 = x2 - x1;
		Vector3s e1 = x0 - x2;
		Vector3s e2 = x0 - x1;
		Vector3s e1_tilde = x3 - x2;
		Vector3s e2_tilde = x3 - x1;
		
		Vector3s n = e0.cross(e1);
		scalar n_length = n.norm();
		Vector3s n_tilde = -e0.cross(e1_tilde);
		scalar n_tilde_length = n_tilde.norm();
		n /= n_length;
		n_tilde /= n_tilde_length;
		
		scalar ka = m_bending_stiffness*3.*e0_rest_sqnorm/restareas;
		
		scalar kb = m_viscous_stiffness*3.*e0_rest_sqnorm/restareas;
		
		scalar rest_phi = m_per_edge_rest_phi(e);
		
		scalar start_phi = m_per_edge_start_phi(e);
		
		scalar dPsi_dTheta;
		
        switch (m_bending_mode) {
            case 0:
                l_compute_dPsi_tantheta(n, n_tilde, e0, rest_phi, ka, start_phi, kb, dPsi_dTheta);
                break;
            case 1:
                l_compute_dPsi_sintheta(n, n_tilde, e0, rest_phi, ka, start_phi, kb, dPsi_dTheta);
                break;
            default:
                l_compute_dPsi_theta(n, n_tilde, e0, rest_phi, ka, start_phi, kb, dPsi_dTheta);
                break;
        }
        
		scalar e0_length = e0.norm();
		
		gradE.segment<3>(idx[0] * 4) += psi_coeff*dPsi_dTheta*(-e0_length/n_length*n.transpose());
		gradE.segment<3>(idx[1] * 4) += psi_coeff*dPsi_dTheta*((-e0/e0_length).dot(e1)/n_length*n.transpose() + (-e0/e0_length).dot(e1_tilde)/n_tilde_length*n_tilde.transpose());
		gradE.segment<3>(idx[2] * 4) += psi_coeff*dPsi_dTheta*((e0/e0_length).dot(e2)/n_length*n.transpose() + (e0/e0_length).dot(e2_tilde)/n_tilde_length*n_tilde.transpose());
		gradE.segment<3>(idx[3] * 4) += psi_coeff*dPsi_dTheta*(-e0_length/n_tilde_length*n_tilde.transpose());
	}
	
	assert(!std::isnan(gradE.sum()));
}

void ShellBendingForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	auto l_bending_stencil_gradientPart = [&] (int e, int idx[4], MatrixXs & dfdx_bending_gradpart, scalar & dPsi_dTheta)
	{
		auto l_compute_dPsi_tantheta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, const scalar rest_phi, const scalar ka, const scalar start_phi, const scalar kb, scalar & dPsi_dTheta, scalar & dPsi_dTheta_dTheta)
		{
			scalar tan_half_theta = (n - n_tilde).norm()/(n + n_tilde).norm();
			scalar sign_angle = n.cross(n_tilde).dot(e0);
			if (std::isnan(sign_angle)) sign_angle = 1.;
			else sign_angle = sign_angle > 0 ? 1. : -1.;
			scalar phi_theta = 2.*sign_angle*tan_half_theta;
			scalar dPhi_dTheta = 1 + tan_half_theta*tan_half_theta;
			dPsi_dTheta = 2.*(ka*(phi_theta - rest_phi)+kb*(phi_theta - start_phi))*dPhi_dTheta;
			dPsi_dTheta_dTheta = 2.*ka*(dPhi_dTheta*dPhi_dTheta + (phi_theta-rest_phi)*tan_half_theta*dPhi_dTheta) +
			2.*kb*(dPhi_dTheta*dPhi_dTheta + (phi_theta-start_phi)*tan_half_theta*dPhi_dTheta);
		};
		
		auto l_compute_dPsi_sintheta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, const scalar rest_phi, const scalar ka, const scalar start_phi, const scalar kb, scalar & dPsi_dTheta, scalar & dPsi_dTheta_dTheta)
		{
			scalar theta = std::atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));
			scalar phi_theta = std::sin(theta/2.);
			dPsi_dTheta = 2.*(ka*(phi_theta - rest_phi)+kb*(phi_theta - start_phi))*(0.5*std::cos(theta/2.));
			dPsi_dTheta_dTheta = 2.*(ka*( (0.5*std::cos(theta/2.))*(0.5*std::cos(theta/2.)) + (phi_theta - rest_phi)*(0.5*0.5*-std::sin(theta/2.)) )+kb*( (0.5*std::cos(theta/2.))*(0.5*std::cos(theta/2.)) + (phi_theta - start_phi)*(0.5*0.5*-std::sin(theta/2.)) ));
		};
		
		auto l_compute_dPsi_theta = [this] (const Vector3s & n, const Vector3s & n_tilde, const Vector3s & e0, const scalar rest_phi, const scalar ka, const scalar start_phi, const scalar kb, scalar & dPsi_dTheta, scalar & dPsi_dTheta_dTheta)
		{
			scalar theta = atan2((n.cross(n_tilde)).dot(e0.normalized()), n.dot(n_tilde));
			scalar phi_theta = theta;
			dPsi_dTheta = 2.*(ka*(phi_theta - rest_phi)+kb*(phi_theta - start_phi));
			dPsi_dTheta_dTheta = 2.*(ka+kb);
		};
		
		const Vector3s& x0 = x.segment<3>(idx[0] * 4);
		const Vector3s& x1 = x.segment<3>(idx[1] * 4);
		const Vector3s& x2 = x.segment<3>(idx[2] * 4);
		const Vector3s& x3 = x.segment<3>(idx[3] * 4);
		
		scalar restareas = m_triangle_rest_area(m_per_unique_edge_triangles(e,0)) + m_triangle_rest_area(m_per_unique_edge_triangles(e,1));
		scalar e0_rest_sqnorm = (m_rest_pos.segment<3>(idx[2] * 4) - m_rest_pos.segment<3>(idx[1] * 4)).squaredNorm();
		
		Vector3s e0 = x2 - x1;
		Vector3s e1 = x0 - x2;
		Vector3s e2 = x0 - x1;
		Vector3s e1_tilde = x3 - x2;
		Vector3s e2_tilde = x3 - x1;
		
		Vector3s n = e0.cross(e1);
		scalar n_length = n.norm();
		Vector3s n_tilde = -e0.cross(e1_tilde);
		scalar n_tilde_length = n_tilde.norm();
		n /= n_length;
		n_tilde /= n_tilde_length;
		
		scalar ka = m_bending_stiffness*3.*e0_rest_sqnorm/restareas;
		
		scalar kb = m_viscous_stiffness*3.*e0_rest_sqnorm/restareas;
		
		scalar rest_phi = m_per_edge_rest_phi(e);
		
		scalar start_phi = m_per_edge_start_phi(e);
		
		scalar dPsi_dTheta_dTheta;
		
        switch (m_bending_mode) {
            case 0:
                l_compute_dPsi_tantheta(n, n_tilde, e0, rest_phi, ka, start_phi, kb, dPsi_dTheta, dPsi_dTheta_dTheta);
                break;
            case 1:
                l_compute_dPsi_sintheta(n, n_tilde, e0, rest_phi, ka, start_phi, kb, dPsi_dTheta, dPsi_dTheta_dTheta);
                break;
            default:
                l_compute_dPsi_theta(n, n_tilde, e0, rest_phi, ka, start_phi, kb, dPsi_dTheta, dPsi_dTheta_dTheta);
                break;
        }
		
		scalar e0_length = e0.norm();
		
		VectorXs dthetadx(12);
		dthetadx.block(0,0,3,1) = -e0_length/n_length*n;
		dthetadx.block(3,0,3,1) = (-e0/e0_length).dot(e1)/n_length*n + (-e0/e0_length).dot(e1_tilde)/n_tilde_length*n_tilde;
		dthetadx.block(6,0,3,1) = (e0/e0_length).dot(e2)/n_length*n + (e0/e0_length).dot(e2_tilde)/n_tilde_length*n_tilde;
		dthetadx.block(9,0,3,1) = -e0_length/n_tilde_length*n_tilde;
		
		dfdx_bending_gradpart = dPsi_dTheta_dTheta*dthetadx*dthetadx.transpose();
	};
	
	auto l_bending_stencil_hessianPart = [&] (int f, int idx[3], MatrixXs & dfdx_bending_hesspart, const VectorXs & e_length, const VectorXs & dPsi_dTheta)
	{
		const Vector3s& x0 = x.segment<3>(idx[0] * 4);
		const Vector3s& x1 = x.segment<3>(idx[1] * 4);
		const Vector3s& x2 = x.segment<3>(idx[2] * 4);
		
		scalar e0_length = e_length(m_per_triangles_unique_edges(f,0));
		scalar e1_length = e_length(m_per_triangles_unique_edges(f,1));
		scalar e2_length = e_length(m_per_triangles_unique_edges(f,2));
		
		Vector3s e0_normalized = (x2-x1)/e0_length;
		Vector3s e1_normalized = (x2-x0)/e1_length;
		Vector3s e2_normalized = (x1-x0)/e2_length;
		
		scalar cos_alpha0 = e1_normalized.dot(e2_normalized);
		scalar cos_alpha1 = -e2_normalized.dot(e0_normalized);
		scalar cos_alpha2 = -e1_normalized.dot(-e0_normalized);
		
		scalar A = 0.5*(x2-x0).cross(x1-x0).norm();
		
		scalar inv_h0 = e0_length/(2.*A);
		scalar inv_h1 = e1_length/(2.*A);
		scalar inv_h2 = e2_length/(2.*A);
		
		Vector3s n = e2_normalized.cross(e1_normalized).normalized();
		
		Vector3s m0 = e0_normalized.cross(n).normalized();
		Vector3s m1 = -e1_normalized.cross(n).normalized();
		Vector3s m2 = e2_normalized.cross(n).normalized();
		
		scalar w[3][3];
		w[0][0] = inv_h0*inv_h0;
		w[0][1] = inv_h0*inv_h1;
		w[0][2] = inv_h0*inv_h2;
		w[1][0] = inv_h1*inv_h0;
		w[1][1] = inv_h1*inv_h1;
		w[1][2] = inv_h1*inv_h2;
		w[2][0] = inv_h2*inv_h0;
		w[2][1] = inv_h2*inv_h1;
		w[2][2] = inv_h2*inv_h2;
		
		Matrix3s M[3];
		M[0] = n*m0.transpose();
		M[1] = n*m1.transpose();
		M[2] = n*m2.transpose();
		
		Matrix3s N0 = M[0]/(e0_length*e0_length);
		Matrix3s N1 = M[1]/(e1_length*e1_length);
		Matrix3s N2 = M[2]/(e2_length*e2_length);
		
		scalar c0 = dPsi_dTheta(m_per_triangles_unique_edges(f,0));
		scalar c1 = dPsi_dTheta(m_per_triangles_unique_edges(f,1));
		scalar c2 = dPsi_dTheta(m_per_triangles_unique_edges(f,2));
		
		scalar d[3];
		d[0] = c2*cos_alpha1 + c1*cos_alpha2 - c0;
		d[1] = c0*cos_alpha2 + c2*cos_alpha0 - c1;
		d[2] = c1*cos_alpha0 + c0*cos_alpha1 - c2;
		
		Matrix3s R[3];
		R[0] = c0*N0;
		R[1] = c1*N1;
		R[2] = c2*N2;
		
		dfdx_bending_hesspart.resize(9,9);
		
		for (int i = 0; i < 3; i++)
		{
			for (int k = i; k < i + 2; k++)
			{
				int j = k % 3;
				
				Matrix3s hij = w[i][j] * (d[i] * M[j].transpose() + d[j] * M[i]);
				if (i == j)
				{
					hij += -R[(i + 1) % 3] - R[(i + 2) % 3];
				}
				else
				{
					// check if we need to transpose. Needs transpose for one of the two triangles of the hinge (doesn't matter which as long as it is consistent throughout sim).
					int e = m_per_triangles_unique_edges(f, (i + 2) % 3); // get edge id
					if (m_per_unique_edge_triangles(e,0) == f)
						hij += R[(i + 2) % 3].transpose();
					else
						hij += R[(i + 2) % 3];
				}
				
				dfdx_bending_hesspart.block(i * 3, j * 3, 3, 3) = hij;
			}
		}
		
		dfdx_bending_hesspart.block(3, 0, 3, 3) = dfdx_bending_hesspart.block(0, 3, 3, 3).transpose();
		dfdx_bending_hesspart.block(0, 6, 3, 3) = dfdx_bending_hesspart.block(6, 0, 3, 3).transpose();
		dfdx_bending_hesspart.block(6, 3, 3, 3) = dfdx_bending_hesspart.block(3, 6, 3, 3).transpose();
	};
	
	VectorXs e_length(m_E_unique.rows());
	VectorXs dPsi_dTheta(m_E_unique.rows());
	dPsi_dTheta.setZero(); // same effect as sigma (boundary boolean)
	
	int base_idx = hessE_index;
	
	threadutils::for_each(0, (int) m_E_unique.rows(), [&] (int e) {
		e_length(e) = (x.segment<3>(m_E_unique(e, 0) * 4) - x.segment<3>(m_E_unique(e, 1) * 4)).norm();
	});
	
	threadutils::for_each(0, (int) m_unique_edge_usable.size(), [&] (int k) {
		const int e = m_unique_edge_usable[k];
		
		assert ( (m_per_unique_edge_triangles(e, 0) != -1) && (m_per_unique_edge_triangles(e, 1) != -1) );
		
		int idx[4];
		idx[0] = m_F(m_per_unique_edge_triangles(e, 0), m_per_unique_edge_triangles_local_corners(e, 0));
		idx[1] = m_E_unique(e, 0);
		idx[2] = m_E_unique(e, 1);
		idx[3] = m_F(m_per_unique_edge_triangles(e, 1), m_per_unique_edge_triangles_local_corners(e, 1));
		const scalar psi_coeff = pow((psi(idx[0]) + psi(idx[1]) + psi(idx[2]) + psi(idx[3])) * 0.25, lambda);
		
		MatrixXs dfdx_bending_gradpart;
		l_bending_stencil_gradientPart(e, idx, dfdx_bending_gradpart, dPsi_dTheta(e));
		
		for(int j = 0; j < 4; ++j) for(int i = 0; i < 4; ++i) for(int s = 0; s < 3; ++s) for(int r = 0; r < 3; ++r)
		{
			int hess_idx = base_idx + k * 16 * 9 + (j * 4 + i) * 9 + s * 3 + r;
			hessE[hess_idx] = Triplets( idx[i] * 4 + r, idx[j] * 4 + s, psi_coeff * dfdx_bending_gradpart(i * 3 + r, j * 3 + s) );
			assert(!std::isnan(hessE[hess_idx].value()));
		}
	});
	
#ifdef USE_FULL_HESSIAN
	base_idx += m_unique_edge_usable.size() * 16 * 9;
	
	threadutils::for_each(0, (int) m_F.rows(), [&] (int f) {
		int idx[3];
		idx[0] = m_F(f,0);
		idx[1] = m_F(f,1);
		idx[2] = m_F(f,2);
		const scalar psi_coeff = pow((psi(idx[0]) + psi(idx[1]) + psi(idx[2])) / 3.0, lambda);
		
		MatrixXs dfdx_bending_hesspart;
		l_bending_stencil_hessianPart(f, idx, dfdx_bending_hesspart, e_length, dPsi_dTheta);
		
		for(int j = 0; j < 3; ++j) for(int i = 0; i < 3; ++i) for(int s = 0; s < 3; ++s) for(int r = 0; r < 3; ++r)
		{
			int hess_idx = base_idx + f * 9 * 9 + (j * 3 + i) * 9 + s * 3 * r;
			hessE[hess_idx] = Triplets( idx[i] * 4 + r, idx[j] * 4 + s, psi_coeff * dfdx_bending_hesspart(i * 3 + r, j * 3 + s) );
			assert(!std::isnan(hessE[hess_idx].value()));
		}
	});
#endif
}

int ShellBendingForce::numHessX()
{
	return m_unique_edge_usable.size() * 16 * 9
#ifdef USE_FULL_HESSIAN
	+ m_F.rows() * 9 * 9
#endif
	;
}

void ShellBendingForce::preCompute()
{
	
}

void ShellBendingForce::updateStartState()
{
	m_per_edge_start_phi.resize(m_per_edge_rest_phi.size());
	computeBendingRestPhi(m_pos, m_per_edge_start_phi);
}

Force* ShellBendingForce::createNewCopy()
{
	return new ShellBendingForce(*this);
}

int ShellBendingForce::flag() const
{
	return 1;
}
