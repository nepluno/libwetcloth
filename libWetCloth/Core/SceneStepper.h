//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
