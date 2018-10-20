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

#ifndef __SCENE_STEPPER__
#define __SCENE_STEPPER__

#include "TwoDScene.h"
#include <functional>
#include <stack>

#include "MathDefs.h"

class SceneStepper
{
public:
	virtual ~SceneStepper();
	
	virtual bool stepScene( TwoDScene& scene, scalar dt ) = 0;
	
	virtual bool stepVelocity( TwoDScene& scene, scalar dt ) = 0;

	virtual bool stepVelocityLagrangian( TwoDScene& scene, scalar dt ) = 0;
    
    virtual bool projectFine( TwoDScene& scene, scalar dt ) = 0;
    
    virtual bool applyPressureDragElasto( TwoDScene& scene, scalar dt ) = 0;
    
    virtual bool applyPressureDragFluid( TwoDScene& scene, scalar dt ) = 0;
	
	virtual bool acceptVelocity( TwoDScene& scene ) = 0;

	virtual bool manifoldPropagate( TwoDScene& scene, scalar dt ) = 0;
    
    virtual bool advectSurfTension( TwoDScene& scene, scalar dt ) = 0;
    
    virtual bool stepImplicitElasto( TwoDScene& scene, scalar dt ) = 0;

    virtual bool stepImplicitElastoLagrangian( TwoDScene& scene, scalar dt ) = 0;
    
    virtual scalar computeDivergence( TwoDScene& scene ) = 0;
    
    virtual void pushFluidVelocity() = 0;
    
    virtual void popFluidVelocity() = 0;
    
    virtual void pushElastoVelocity() = 0;
    
    virtual void popElastoVelocity() = 0;
	
	virtual bool advectScene( TwoDScene& scene, scalar dt );
	
	virtual std::string getName() const = 0;
	
	virtual void setUseApic(bool apic);

	virtual bool useApic() const;
	
	// tools function
	void mapNodeToSoftParticles( const TwoDScene& scene, const std::vector< VectorXs >& node_vec_x, const std::vector< VectorXs >& node_vec_y, const std::vector< VectorXs >& node_vec_z, VectorXs& part_vec ) const;
    
    // tools function
    void buildNodeToSoftParticlesMat( const TwoDScene& scene,
                                     const std::vector< VectorXi >& node_global_indices_x,
                                     const std::vector< VectorXi >& node_global_indices_y,
                                     const std::vector< VectorXi >& node_global_indices_z,
                                     const std::vector< Vector3i >& effective_node_indices,
                                     TripletXs& tri_W,
                                     SparseXs& W ) const;
    
    void buildLocalGlobalMapping( const TwoDScene& scene,
                                 std::vector< VectorXi >& node_global_indices_x,
                                 std::vector< VectorXi >& node_global_indices_y,
                                 std::vector< VectorXi >& node_global_indices_z,
                                 std::vector< Vector3i >& effective_node_indices,
                                 std::vector< Vector3i >& dof_ijk);
	
	void mapSoftParticlesToNode( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z, const VectorXs& part_vec ) const;
	
	void mapSoftParticlesToNodeSqr( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z, const VectorXs& part_vec ) const;

	void allocateNodeVectors( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z) const;
    
    void allocateNodeVectors( const TwoDScene& scene, std::vector< VectorXi >& node_vec_x, std::vector< VectorXi >& node_vec_y, std::vector< VectorXi >& node_vec_z) const;
    
    void allocateCenterNodeVectors( const TwoDScene& scene, std::vector< VectorXi >& node_vec_p ) const;
    
    void allocateCenterNodeVectors( const TwoDScene& scene, std::vector< VectorXs >& node_vec_p ) const;
	
	scalar dotNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az, const std::vector< VectorXs >& node_vec_bx, const std::vector< VectorXs >& node_vec_by, const std::vector< VectorXs >& node_vec_bz ) const;

	scalar dotNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az, const std::vector< VectorXs >& node_vec_bx, const std::vector< VectorXs >& node_vec_by, const std::vector< VectorXs >& node_vec_bz, const VectorXs& twist_vec_a, const VectorXs& twist_vec_b ) const;
	
	scalar dotNodeVectors( const std::vector< VectorXs >& node_vec_a, const std::vector< VectorXs >& node_vec_b ) const;
	
	scalar lengthNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az, const VectorXs& twist_vec ) const;

	scalar lengthNodeVectors( const std::vector< VectorXs >& node_vec_ax, const std::vector< VectorXs >& node_vec_ay, const std::vector< VectorXs >& node_vec_az ) const;
	
	scalar lengthNodeVectors( const std::vector< VectorXs >& node_vec ) const;
	
	void mapGaussToNode( const TwoDScene& scene, std::vector< VectorXs >& node_vec_x, std::vector< VectorXs >& node_vec_y, std::vector< VectorXs >& node_vec_z, const MatrixXs& gauss_vec ) const;
    
    void allocateLagrangianVectors( const TwoDScene& scene, VectorXs& vec );
protected:
	bool m_apic;
};

#endif
