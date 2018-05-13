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

#ifndef __STRAND_FORCE_H__
#define __STRAND_FORCE_H__

#include <Eigen/Core>
#include <iostream>
#include <unordered_set>
#include "../Force.h"

#include "Definitions.h"
#include "../TwoDScene.h"
#include "Dependencies/ElasticityParameters.h"
#include "Dependencies/DegreesOfFreedom.h"
#include "Dependencies/ReferenceFrames.h"
#include "Dependencies/MaterialFrames.h"
#include "Dependencies/Kappas.h"
#include "Dependencies/Twists.h"
#include "Dependencies/BendingProducts.h"
#include "StrandParameters.h"

struct StrandState
{
	StrandState( const VecX& initDofs, BendingMatrixBase& bendingMatrixBase );
	
	DOFs m_dofs;
	Edges m_edges;
	Lengths m_lengths;
	Tangents m_tangents;
	mutable ReferenceFrames1 m_referenceFrames1;
	mutable ReferenceFrames2 m_referenceFrames2;
	ReferenceTwists m_referenceTwists;
	Twists m_twists;
	CurvatureBinormals m_curvatureBinormals;
	TrigThetas m_trigThetas;
	mutable MaterialFrames<1> m_materialFrames1;
	mutable MaterialFrames<2> m_materialFrames2;
	Kappas m_kappas;
	GradKappas m_gradKappas;
	GradTwists m_gradTwists;
	GradTwistsSquared m_gradTwistsSquared;
	HessKappas m_hessKappas;
	HessTwists m_hessTwists;
	BendingProducts m_bendingProducts;
};

struct StartState
{ // used for Viscous updates that depend of start of step state
	StartState( const VecX& initDofs );
	
	DOFs m_dofs;
	Edges m_edges;
	Lengths m_lengths;
	Tangents m_tangents;
	mutable ReferenceFrames1 m_referenceFrames1;
	mutable ReferenceFrames2 m_referenceFrames2;
	ReferenceTwists m_referenceTwists;
	Twists m_twists;
	CurvatureBinormals m_curvatureBinormals;
	TrigThetas m_trigThetas;
	mutable MaterialFrames<1> m_materialFrames1;
	mutable MaterialFrames<2> m_materialFrames2;
	Kappas m_kappas;
};

class StrandForce : public Force
{
public:
	
	StrandForce( const std::shared_ptr<TwoDScene>& scene, const std::vector<int>& consecutiveVertices, const int& parameterIndex, int globalIndex );
	
	virtual ~StrandForce();
	
	virtual Force* createNewCopy();
	
	virtual void preCompute();
	
	virtual void updateStartState();
	
	virtual void addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E );
	
	virtual void addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE );
	
	virtual void addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt );
	
	virtual int numHessX();
	
	virtual int flag() const;
	
	virtual const char* name();
	
	int getGlobalIndex() const { return m_globalIndex; }
	int getNumVertices() const { return (int) m_verts.size(); }
	
	Vec2Array& alterRestKappas()
	{
		return m_restKappas;
	}
	

    void updateStrandState();
	
	// private: //todo need to make get/set calls for twodscene in order to make this private again.
	
	int getNumEdges() const { return (int) m_verts.size() - 1; }
	
	void resizeInternals();
	void freezeRestShape( unsigned begin, unsigned end, scalar damping = 0. );
	void updateRestShape( const VecX& x_restshape, scalar damping = 0. );
	
	void updateEverythingThatDependsOnRestLengths();
	
	int numConstraintNonViscous();
	int numConstraintViscous();
	
	void recomputeGlobal();
	void clearStored();
	void getLocalAffectedVars( int colidx, std::vector< std::pair<int,int> >& vars );
	
	template<typename AccumulatedT>
	void accumulateQuantity( AccumulatedT& accumulated );
	
	//// FOSSSim related //////////////////////////////////////////////////
	std::vector< int > m_verts; // in order root to tip
	VectorXs m_packing_fraction;
	int m_globalIndex; // Global index amongst the hairs
	std::shared_ptr<StrandParameters> m_strandParams;
	std::shared_ptr<TwoDScene> m_scene;
	bool m_requiresExactForceJacobian;
	
	// increase memory, reduce re-computation
	scalar m_strandEnergyUpdate;
	VecX m_strandForceUpdate;
	TripletXs m_strandHessianUpdate;
	
	//// Strand State (implicitly the end of timestep state, evolved from rest config) ////////////////////////
	StrandState* m_strandState; // future state
	StartState* m_startState; // current state
	
	//// Rest shape //////////////////////////////////////////////////////
	std::vector<scalar> m_restLengths; // The following four members depend on m_restLengths, which is why updateEverythingThatDependsOnRestLengths() must be called
	scalar m_totalRestLength;
	std::vector<scalar> m_VoronoiLengths; // rest length around each vertex
	std::vector<scalar> m_invVoronoiLengths; // their inverses
	std::vector<scalar> m_vertexMasses;
	Vec2Array m_restKappas;
	std::vector<scalar> m_restTwists;
	
	//// Friends /////////////////////////////////////////////////////////
	friend class Viscous;
	friend class NonViscous;
	template<typename ViscousT> friend class StretchingForce;
	template<typename ViscousT> friend class BendingForce;
	template<typename ViscousT> friend class TwistingForce;
};

#endif
