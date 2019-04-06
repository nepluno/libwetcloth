//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2016 Gabriel Cirio
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

#ifndef __SHELL_BENDING_FORCE_H__
#define __SHELL_BENDING_FORCE_H__

#include <Eigen/Core>

#include "../../MathDefs.h"
#include "../../Force.h"

class ShellBendingForce : public Force
{
protected:
	const VectorXs & m_pos;
	const VectorXs & m_rest_pos;
	const MatrixXi & m_F;
	const VectorXs & m_triangle_rest_area;
	
	const MatrixXi & m_E_unique;
	
	const MatrixXi & m_per_unique_edge_triangles;
	const MatrixXi & m_per_unique_edge_triangles_local_corners;
	const MatrixXi & m_per_triangles_unique_edges;
	
	const scalar& m_young_modulus;
	const scalar& m_viscous_modulus;
	const scalar& m_poisson_ratio;
	
	scalar m_thickness;
	
	bool m_apply_viscous;
	
	scalar m_bending_stiffness;
	scalar m_viscous_stiffness;
	VectorXs m_per_edge_rest_phi; // theta or tantheta, depending on choice of bending formulation
	VectorXs m_per_edge_start_phi;
	std::vector<int> m_unique_edge_usable;
	VectorXs m_multipliers;
	
	void computeBendingRestPhi(const VectorXs & rest_pos, VectorXs & rest_phis);
    
    int m_bending_mode;
	
public:
	
	ShellBendingForce(const VectorXs & pos,
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
                      int bending_mode,
					  bool apply_viscous);
	
	virtual ~ShellBendingForce();
	
	virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E );
	
	virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE );
	
	virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt );

	virtual void updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt );
	
	virtual int numHessX();
	
	virtual void preCompute();
	
	virtual void updateStartState();
	
	virtual Force* createNewCopy();
	
	virtual int flag() const;
};

#endif
