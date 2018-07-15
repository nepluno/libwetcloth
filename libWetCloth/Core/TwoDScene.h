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

#ifndef __TWO_D_SCENE_H__
#define __TWO_D_SCENE_H__

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <fstream>
#include "Force.h"
#include "DER/StrandParameters.h"
#include "sorter.h"
#include "Script.h"
#include "DistanceFields.h"

class StrandForce;
class AttachForce;

#ifdef PC_NONE
#undef PC_NONE
#endif

enum ParticleClassifier
{
    PC_NONE,
    PC_S,
    PC_s,
    PC_o,
    PC_l,
    PC_L
};

enum NODE_STATE
{
	NS_NONE,
	NS_FLUID,
	NS_SOLID
};

struct LiquidInfo
{
	scalar liquid_density;
	scalar air_density;
	scalar surf_tension_coeff;
	scalar viscosity;
	scalar air_viscosity;
	scalar rest_contact_angle;
	scalar yazdchi_power;
	scalar pore_radius;
    scalar yarn_diameter;
    scalar rest_volume_fraction;
	scalar lambda;
    scalar cohesion_coeff;
    scalar correction_multiplier;
    scalar correction_strength;
    scalar flip_coeff;
    scalar elasto_flip_asym_coeff;
    scalar elasto_flip_coeff;
    scalar elasto_advect_coeff;
	scalar particle_cell_multiplier;
    scalar levelset_young_modulus;
	scalar liquid_boundary_friction;
	scalar levelset_thickness;
	scalar elasto_capture_rate;
    int correction_step;
    int bending_scheme;
    int iteration_print_step;
    bool use_cohesion;
	bool solid_cohesion;
	bool soft_cohesion;
    bool solve_solid;
    bool use_nonlinear_drag;
    bool use_drag;
    bool apply_pressure_solid;
    bool use_levelset_force;
    bool apply_pressure_manifold;
    bool use_twist;
    bool use_bicgstab;
    bool use_amgpcg_solid;
    bool use_pcr;
    bool apply_pore_pressure_solid;
    bool propagate_solid_velocity;
    bool check_divergence;
	bool use_varying_fraction;
    bool compute_viscosity;
    bool implicit_viscosity;
    bool drag_by_future_solid;
    bool drag_by_air;
    bool init_nonuniform_fraction;
    bool use_group_precondition;
    bool use_lagrangian_mpm;
    bool use_cosolve_angular;
	
    friend std::ostream& operator<<(std::ostream&, const LiquidInfo&);
};

std::ostream& operator<<(std::ostream&, const LiquidInfo&);

struct RayTriInfo
{
    int start_geo_id;
    Vector3s end;
    Vector3s norm;
    int intersect_geo_id;
    Vector2s uv;
    scalar dist;
    scalar volume_frac;
    scalar c0;
    scalar c1;
	scalar weight;
};

class TwoDScene : public std::enable_shared_from_this<TwoDScene>
{
    const static int m_kernel_order = 2;
    const static int m_num_armor = 0;
    
public:
    
    TwoDScene();
    
    TwoDScene( const TwoDScene& otherscene ) = delete;
    
    ~TwoDScene();
    
    int getNumParticles() const;
    int getNumEdges() const;
	int getNumFaces() const;
	int getNumSurfels() const;
	int getNumGausses() const;
	int getNumBuckets() const;

	int getNumStrandParameters() const;
    
    const std::vector<int> getParticleGroup() const;
    
    const std::vector<unsigned char>& getFixed() const;
    
    const VectorXs& getX() const;
    
    VectorXs& getX();
    
    const VectorXs& getV() const;
    
    VectorXs& getV();
	
	const VectorXs& getFluidV() const;
	
	VectorXs& getFluidV();
	
    const VectorXs& getM() const;
    
    VectorXs& getM();
	
	const VectorXs& getFluidM() const;
	
	VectorXs& getFluidM();
	
	const VectorXs& getVol() const;
	
	VectorXs& getVol();
	
	const VectorXs& getFluidVol() const;
	
	VectorXs& getFluidVol();
    
    LiquidInfo& getLiquidInfo();
	
	const LiquidInfo& getLiquidInfo() const;
	
	const std::vector< VectorXs >& getParticleDiv() const;
	
	const std::vector< int >& getParticleEdges(int pidx) const;
	
	const std::vector< std::pair<int, scalar> >& getParticleFaces(int pidx) const;
	
	const VectorXs& getVolumeFraction() const;
	
	VectorXs& getVolumeFraction();
	
	const VectorXs& getRadius() const;
	
	VectorXs& getRadius();
	
	const VectorXs& getGaussX() const;
	
	VectorXs& getGaussX();
	
	const VectorXs& getGaussV() const;
	
	VectorXs& getGaussV();
    
    const MatrixXs& getGaussNormal() const;
    
    MatrixXs& getGaussNormal();
	
	const VectorXs& getGaussDV() const;
	
	VectorXs& getGaussDV();
	
	const VectorXs& getGaussFluidM() const;
	
	VectorXs& getGaussFluidM();
	
	const VectorXs& getGaussFluidVol() const;
	
	VectorXs& getGaussFluidVol();
	
	const VectorXs& getGaussVolumeFraction() const;
	
	VectorXs& getGaussVolumeFraction();
	
	const VectorXs& getGaussM() const;
	
	VectorXs& getGaussM();
	
	const VectorXs& getGaussVol() const;
	
	VectorXs& getGaussVol();
	
	const MatrixXs& getGaussFe() const;
	
	MatrixXs& getGaussFe();
	
	const MatrixXs& getGaussd() const;
	
	MatrixXs& getGaussd();
    
    const MatrixXs& getGaussDinv() const;
    
    MatrixXs& getGaussDinv();
	
    const MatrixXs& getGaussD() const;
    
    MatrixXs& getGaussD();
    
    const Matrix27x3s& getGaussdFdXX(int pidx) const;

    const Matrix27x3s& getGaussdFdXY(int pidx) const;
	
    const Matrix27x3s& getGaussdFdXZ(int pidx) const;
	
	int getDefaultNumNodes() const;
    
	int getNumNodesX(int bucket_idx) const;
	
	int getNumNodesY(int bucket_idx) const;
	
	int getNumNodesZ(int bucket_idx) const;
	
	int getNumNodesP(int bucket_idx) const;
	
	int getNumNodesSolidPhi(int bucket_idx) const;
	
	inline bool isFluid(int pidx) const;
	
	const std::vector< std::shared_ptr< DistanceField > >& getGroupDistanceField() const;
	
	const std::vector< VectorXuc >& getNodeLiquidValidX() const;
	
	std::vector< VectorXuc >& getNodeLiquidValidX();
	
	const std::vector< VectorXuc >& getNodeLiquidValidY() const;
	
	std::vector< VectorXuc >& getNodeLiquidValidY();
	
	const std::vector< VectorXuc >& getNodeLiquidValidZ() const;
	
	std::vector< VectorXuc >& getNodeLiquidValidZ();
	
	const std::vector< VectorXuc >& getNodeStateX() const;
	
	std::vector< VectorXuc >& getNodeStateX();
	
	const std::vector< VectorXuc >& getNodeStateY() const;
	
	std::vector< VectorXuc >& getNodeStateY();
	
	const std::vector< VectorXuc >& getNodeStateZ() const;
	
	std::vector< VectorXuc >& getNodeStateZ();
	
	const std::vector< VectorXi >& getNodeCompressedIndexX() const;
	
	const std::vector< VectorXi >& getNodeCompressedIndexY() const;
	
	const std::vector< VectorXi >& getNodeCompressedIndexZ() const;
    
    const std::vector< VectorXi >& getNodeCompressedIndexP() const;
    
	const VectorXs& getNodePosX(int bucket_idx) const;
	
	VectorXs& getNodePosX(int bucket_idx);
	
	const VectorXs& getNodePosY(int bucket_idx) const;
	
	VectorXs& getNodePosY(int bucket_idx);
	
	const VectorXs& getNodePosZ(int bucket_idx) const;
	
	VectorXs& getNodePosZ(int bucket_idx);
    
    const VectorXs& getNodePosEX(int bucket_idx) const;
    
    VectorXs& getNodePosEX(int bucket_idx);
    
    const VectorXs& getNodePosEY(int bucket_idx) const;
    
    VectorXs& getNodePosEY(int bucket_idx);
    
    const VectorXs& getNodePosEZ(int bucket_idx) const;
    
    VectorXs& getNodePosEZ(int bucket_idx);
	
	const VectorXs& getNodePosP(int bucket_idx) const;
	
	VectorXs& getNodePosP(int bucket_idx);
    
    const VectorXi& getNodeIndicesP(int bucket_idx) const;
    
    VectorXi& getNodeIndicesP(int bucket_idx);
	
	const VectorXs& getNodePosSolidPhi(int bucket_idx) const;
	
	VectorXs& getNodePosSolidPhi(int bucket_idx);
    
    const std::vector<VectorXs>& getNodePosSolidPhi() const;
    
    std::vector<VectorXs>& getNodePosSolidPhi();
    
    const std::vector<VectorXi>& getNodeIndicesX() const;
    
    std::vector<VectorXi>& getNodeIndicesX();
    
    const std::vector<VectorXi>& getNodeIndicesY() const;
    
    std::vector<VectorXi>& getNodeIndicesY();
    
    const std::vector<VectorXi>& getNodeIndicesZ() const;
    
    std::vector<VectorXi>& getNodeIndicesZ();
    
    const std::vector<VectorXi>& getNodeIndicesP() const;
    
    std::vector<VectorXi>& getNodeIndicesP();
    
    const std::vector<VectorXs>& getNodeSolidPhi() const;
    
    std::vector<VectorXs>& getNodeSolidPhi();
	
    const std::vector< VectorXs >& getNodeOrientationX() const;
    
    const std::vector< VectorXs >& getNodeOrientationY() const;
    
    const std::vector< VectorXs >& getNodeOrientationZ() const;
    
    const std::vector< VectorXs >& getNodeShapeFactorX() const;
    
    const std::vector< VectorXs >& getNodeShapeFactorY() const;
    
    const std::vector< VectorXs >& getNodeShapeFactorZ() const;
    
	const std::vector< VectorXs >& getNodePressure() const;
	
	std::vector< VectorXs >& getNodePressure();
	
	const std::vector< VectorXs >& getNodeCellSolidPhi() const;
	
	std::vector< VectorXs >& getNodeCellSolidPhi();
	
	const std::vector< VectorXs >& getNodeVelocityX() const;
	
	std::vector< VectorXs >& getNodeVelocityX();
	
	const std::vector< VectorXs >& getNodeVelocityY() const;
	
	std::vector< VectorXs >& getNodeVelocityY();
	
	const std::vector< VectorXs >& getNodeVelocityZ() const;
	
	std::vector< VectorXs >& getNodeVelocityZ();
	
	const std::vector< VectorXs >& getNodeFluidVelocityX() const;
	
	std::vector< VectorXs >& getNodeFluidVelocityX();
	
	const std::vector< VectorXs >& getNodeFluidVelocityY() const;
	
	std::vector< VectorXs >& getNodeFluidVelocityY();
	
	const std::vector< VectorXs >& getNodeFluidVelocityZ() const;
	
	std::vector< VectorXs >& getNodeFluidVelocityZ();
	
	const std::vector< VectorXs >& getNodePorePressureP() const;
	
	std::vector< VectorXs >& getNodePorePressureP();
	
	const std::vector< VectorXs >& getNodeSaturationP() const;
	
	std::vector< VectorXs >& getNodeSaturationP();
	
	const std::vector< VectorXs >& getNodeSaturationX() const;
	
	std::vector< VectorXs >& getNodeSaturationX();
	
	const std::vector< VectorXs >& getNodeSaturationY() const;
	
	std::vector< VectorXs >& getNodeSaturationY();
	
	const std::vector< VectorXs >& getNodeSaturationZ() const;
	
	std::vector< VectorXs >& getNodeSaturationZ();
	
	const std::vector< VectorXs >& getNodePsiX() const;
	
	std::vector< VectorXs >& getNodePsiX();
	
	const std::vector< VectorXs >& getNodePsiY() const;
	
	std::vector< VectorXs >& getNodePsiY();
	
	const std::vector< VectorXs >& getNodePsiZ() const;
	
	std::vector< VectorXs >& getNodePsiZ();
	
	const std::vector< VectorXs >& getNodeMassX() const;
	
	std::vector< VectorXs >& getNodeMassX();
	
	const std::vector< VectorXs >& getNodeMassY() const;
	
	std::vector< VectorXs >& getNodeMassY();
	
	const std::vector< VectorXs >& getNodeMassZ() const;
	
	std::vector< VectorXs >& getNodeMassZ();
	
	const std::vector< VectorXs >& getNodeFluidMassX() const;
	
	std::vector< VectorXs >& getNodeFluidMassX();
	
	const std::vector< VectorXs >& getNodeFluidMassY() const;
	
	std::vector< VectorXs >& getNodeFluidMassY();
	
	const std::vector< VectorXs >& getNodeFluidMassZ() const;
	
	std::vector< VectorXs >& getNodeFluidMassZ();
	
	const std::vector< VectorXs >& getNodeFluidVolX() const;
	
	std::vector< VectorXs >& getNodeFluidVolX();
	
	const std::vector< VectorXs >& getNodeFluidVolY() const;
	
	std::vector< VectorXs >& getNodeFluidVolY();
	
	const std::vector< VectorXs >& getNodeFluidVolZ() const;
	
	std::vector< VectorXs >& getNodeFluidVolZ();
	
	const std::vector< VectorXs >& getNodeVolX() const;
	
	std::vector< VectorXs >& getNodeVolX();
	
	const std::vector< VectorXs >& getNodeVolY() const;
	
	std::vector< VectorXs >& getNodeVolY();
	
	const std::vector< VectorXs >& getNodeVolZ() const;
	
	std::vector< VectorXs >& getNodeVolZ();
	
	const std::vector< VectorXs >& getNodeSolidWeightX() const;
	
	const std::vector< VectorXs >& getNodeSolidWeightY() const;
	
	const std::vector< VectorXs >& getNodeSolidWeightZ() const;
	
	const std::vector< VectorXs >& getNodeSolidVelX() const;
	
	const std::vector< VectorXs >& getNodeSolidVelY() const;
	
	const std::vector< VectorXs >& getNodeSolidVelZ() const;
	
	const std::vector< VectorXs >& getNodeLiquidPhi() const;
	
	const std::vector< VectorXs >& getNodeCurvatureP() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracCentre() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracU() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracV() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracW() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracEX() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracEY() const;
    
    const std::vector< VectorXs >& getNodeLiquidVolFracEZ() const;
    
    const std::vector< VectorXi >& getNodeIndexEdgeX() const;
    
    const std::vector< VectorXi >& getNodeIndexEdgeY() const;
    
    const std::vector< VectorXi >& getNodeIndexEdgeZ() const;
    
	const Matrix27x3i& getParticleNodesX( int pidx ) const;
	
	Matrix27x3i& getParticleNodesX( int pidx );
	
	const Matrix27x3i& getGaussNodesX( int pidx ) const;
	
	Matrix27x3i& getGaussNodesX( int pidx );
	
	const Matrix27x3i& getParticleNodesY( int pidx ) const;
	
	Matrix27x3i& getParticleNodesY( int pidx );
	
	const Matrix27x3i& getGaussNodesY( int pidx ) const;
	
	Matrix27x3i& getGaussNodesY( int pidx );
	
	const Matrix27x3i& getParticleNodesZ( int pidx ) const;
	
	Matrix27x3i& getParticleNodesZ( int pidx );
    
    const Matrix27x3i& getParticleNodesSolidPhi( int pidx ) const;
    
    Matrix27x3i& getParticleNodesSolidPhi( int pidx );
    
    const Matrix27x3s& getParticleGradsSolidPhi( int pidx ) const;
    
    Matrix27x3s& getParticleGradsSolidPhi( int pidx );
	
	const Matrix27x3i& getGaussNodesZ( int pidx ) const;
	
	Matrix27x3i& getGaussNodesZ( int pidx );

    const std::vector< Matrix27x4s >& getParticleWeights() const;
	
	const Matrix27x4s& getParticleWeights( int pidx ) const;
	
	Matrix27x4s& getParticleWeights( int pidx );
    
    const Matrix27x3s& getGaussGradsX(int pidx) const;
    
    Matrix27x3s& getGaussGradsX(int pidx);
    
    const Matrix27x3s& getParticleGradsX(int pidx) const;
    
    Matrix27x3s& getParticleGradsX(int pidx);
	
	const Matrix27x3s& getGaussGradsY(int pidx) const;
	
	Matrix27x3s& getGaussGradsY(int pidx);
	
	const Matrix27x3s& getParticleGradsY(int pidx) const;
	
	Matrix27x3s& getParticleGradsY(int pidx);
	
	const Matrix27x3s& getGaussGradsZ(int pidx) const;
	
	Matrix27x3s& getGaussGradsZ(int pidx);
	
	const Matrix27x3s& getParticleGradsZ(int pidx) const;
	
	Matrix27x3s& getParticleGradsZ(int pidx);
	
	const std::vector< VectorXi >& getPressureNeighbors() const;
	
	const std::vector< VectorXi >& getNodePressureIndexX() const;
	
	const std::vector< VectorXi >& getNodePressureIndexY() const;
	
	const std::vector< VectorXi >& getNodePressureIndexZ() const;
	
	void distributeFluidElasto(const scalar& dt);
	
	void distributeElastoFluid();
	
	void swapParticles(int i, int j);
	
    void mapParticleNodesAPIC(); // particles to nodes mapping

	void mapParticleSaturationPsiNodes();
	
	void updatePorePressureNodes();
	
    void mapNodeParticlesAPIC(); // nodes to particles mapping

    void resizeParticleSystem( int num_particles );
    
    void conservativeResizeParticles(int num_particles);
	
	void conservativeResizeEdges(int num_edges);
	
	void conservativeResizeFaces(int num_faces);
	
	void updateGaussSystem(scalar dt);
	
	void updateGaussManifoldSystem();
	
	void updateGaussAccel();
	
	void updateManifoldOperators();
	
	void updateSolidVolumeFraction();
	
	void initGaussSystem();
    
    VectorXs getPosition(int particle);
    
    int getDof(int particle) const;
    
    void setPosition( int particle, const Vector3s& pos );
	
	void setLiquidInfo( const LiquidInfo& info );
    
    void setTheta(int particle, const scalar theta);
    
    void setVelocity( int particle, const Vector3s& vel );
    
    void setOmega( int particle, const scalar omega );
    
    void setVolume(int particle, const scalar& volume);
	
    void setFluidVolume(int particle, const scalar& volume);
	
    void setVolumeFraction(int particle, const scalar& vol_frac);
	
	void setGroup(int particle, int group);
	
    void setMass( int particle, const scalar& mass, const scalar& second_moments );
	
    void setFluidMass( int particle, const scalar& mass, const scalar& second_moments );
	
    void setRadius( int particle, const scalar& radiusA, const scalar& radiusB );
	
    void setFixed( int particle, unsigned char fixed );
	
    void setTwist( int particle, bool twist );
    
    void setEdge( int idx, const std::pair<int,int>& edge );
	
	void setFace( int idx, const Vector3i& face );
	
	void setParticleToParameters( int idx, int params );
    
    scalar getParticleRestLength( int idx ) const;
	
	scalar getParticleRestArea( int idx ) const;
	
	scalar getCellSize() const;
	
	scalar getInverseDCoeff() const;

    scalar getGaussDensity(int pidx) const;
	
	scalar getInitialVolumeFraction(int pidx) const;
	
	scalar getCapillaryPressure(const scalar& psi) const;
	
	scalar getGaussRadius(int pidx, int dir) const;
    
    scalar getMu(int pidx) const;
    
    scalar getLa(int pidx) const;

   	scalar getViscousModulus(int pidx) const;
	
	scalar getYoungModulus(int pidx) const;

	scalar getShearModulus(int pidx) const;
    
    scalar getAttachMultiplier(int pidx) const;

	scalar getCollisionMultiplier(int pidx) const;

	void loadAttachForces();
	
	Vector3s getTwistDir( int particle ) const;

	Vector3s getRestTwistDir( int particle ) const;
    
    void setEdgeRestLength( int idx, const scalar& l0 );
    
    void setEdgeToParameter( int idx, int params );
	
	void setFaceRestArea( int idx, const scalar& a0 );
    
    void setFaceToParameter( int idx, int params );
	
	const VectorXs& getFaceRestArea() const;
	
	const VectorXs& getEdgeRestLength() const;
    
    unsigned char isFixed( int particle ) const;
    
    bool isGaussFixed(int pidx) const;
	
    bool isTwist( int particle ) const;
    
    const std::vector<bool>& getTwist() const;

    void clearEdges();
    
    const std::vector<int>& getSurfels() const;
    
    const MatrixXi& getEdges() const;

    const MatrixXi& getFaces() const;
	
	void insertStrandParameters( const std::shared_ptr<StrandParameters>& newparams );

    
    std::shared_ptr<StrandParameters>& getStrandParameters( const int index );
    
    const Vector2iT getEdge(int edg) const;
    
    void insertForce( const std::shared_ptr<Force>& newforce );
    
    void setTipVerts( int particle, bool tipVerts );
    
    void accumulateGradU( VectorXs& F, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs());

    void accumulateFluidGradU( VectorXs& F, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs());
    
    void accumulateFluidNodeGradU( std::vector< VectorXs >& node_rhs_x, std::vector< VectorXs >& node_rhs_y, std::vector< VectorXs >& node_rhs_z, scalar coeff );

    void accumulateManifoldFluidGradU( VectorXs& F );
	
	void accumulateManifoldGradPorePressure( VectorXs& F );
	
    void accumulateGaussGradU(MatrixXs& F, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs());
    
    void precompute();

    void postcompute(VectorXs& v, const scalar& dt);

	void stepScript(const scalar& dt, const scalar& current_time);
	
	void applyScript(const scalar& dt);
	
	void updateStartState();
    
    void accumulateddUdxdx( TripletXs& A, const scalar& dt, int base_idx, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs() );
	
    void accumulateAngularddUdxdx( TripletXs& A, const scalar& dt, int base_idx, const VectorXs& dx = VectorXs(), const VectorXs& dv = VectorXs() );

    void computedEdFe(MatrixXs& dFe);
	
    scalar computeKineticEnergy() const;
    scalar computePotentialEnergy() const;
    scalar computeTotalEnergy() const;
    
    void checkConsistency();

	bool isTip(int particle) const;
	
	int getNumBucketColors() const;
	
	const std::vector<int>& getFluidIndices() const;
	
	std::vector<int>& getFluidIndices();
	
	int getNumFluidParticles() const;
	
	int getNumElastoParticles() const;
    
    int getNumSoftElastoParticles() const;
	
	Sorter& getParticleBuckets();
	
	const Sorter& getParticleBuckets() const;
	
	Sorter& getGaussBuckets();
	
	const Sorter& getGaussBuckets() const;
	
	void updateParticleBoundingBox();
	void rebucketizeParticles();
	void resampleNodes();
	void updateElastoParticleWeights(scalar dt);
	void updateLiquidParticleWeights(scalar dt);
	void updateParticleWeights(scalar dt, int start, int end);
	void updateGaussWeights(scalar dt);
	void postProcess(scalar dt);
	void updateSolidWeights();
	void updateLiquidPhi(scalar dt);
    void estimateVolumeFractions(std::vector< VectorXs >& volumes, const std::vector< VectorXs >& node_pos);
    void updateOptiVolume();
    void splitLiquidParticles();
    void mergeLiquidParticles();
    void relabelLiquidParticles();
	void updateRestPos();
    bool useBiCGSTAB() const;
    bool useAMGPCGSolid() const;
    bool propagateSolidVelocity() const;
	
	void correctLiquidParticles(const scalar& dt);
	
	const VectorXs& getRestPos() const;

	 VectorXs& getRestPos() ;	
	
	void initGroupPos();
    
    void updateDeformationGradient(scalar dt);
    void updatePlasticity(scalar dt);
    void updateGaussdFdx();
	
	void updateTotalMass();

	void setBucketInfo( const scalar& bucket_size, int num_cells, int kernel_order );
	
	void buildNodeParticlePairs();
	
	void insertScript( const std::shared_ptr< Script >& script );
    
    const std::vector< std::shared_ptr<AttachForce> >& getAttachForces() const;
	
	const std::vector< std::pair<int, int> >& getNodeParticlePairsX(int bucket_idx, int pidx) const;

	const std::vector< std::pair<int, int> >& getNodeGaussPairsX(int bucket_idx, int pidx) const;
	
	const std::vector< std::pair<int, int> >& getNodeParticlePairsY(int bucket_idx, int pidx) const;
	
	const std::vector< std::pair<int, int> >& getNodeGaussPairsY(int bucket_idx, int pidx) const;
	
	const std::vector< std::pair<int, int> >& getNodeParticlePairsZ(int bucket_idx, int pidx) const;
	
	const std::vector< std::pair<int, int> >& getNodeGaussPairsZ(int bucket_idx, int pidx) const;
    
    scalar getFrictionAlpha(int pidx) const;
    scalar getFrictionBeta(int pidx) const;
	
	Eigen::Quaternion<scalar>& getGroupRotation(int group_idx);
	Vector3s& getGroupTranslation(int group_idx);
	
	Eigen::Quaternion<scalar>& getPrevGroupRotation(int group_idx);
	Vector3s& getPrevGroupTranslation(int group_idx);
	
	void resizeGroups(int num_group);
	
	std::shared_ptr< DistanceField >& getGroupDistanceField(int igroup);
	
	const std::shared_ptr< DistanceField >& getGroupDistanceField(int igroup) const;
	
	std::vector< std::shared_ptr<DistanceField> >& getDistanceFields();
	
	const std::vector< std::shared_ptr<DistanceField> >& getDistanceFields() const;
    
    const VectorXuc& getOutsideInfo() const;
	
	void sampleSolidDistanceFields();
	
	void sampleLiquidDistanceFields( scalar cur_time );
	
	scalar computePhiVel(const Vector3s& pos, Vector3s& vel, const std::function< bool(const std::shared_ptr<DistanceField>&) > selector = nullptr) const;
	
	scalar computePhi(const Vector3s& pos, const std::function< bool(const std::shared_ptr<DistanceField>&) > selector = nullptr) const;
    
    void dump_geometry(std::string filename);
	
	int getKernelOrder() const;
	
	void updateSolidPhi();

	void updateMultipliers(const scalar& dt);
	
	void solidProjection(const scalar& dt);
    
    void terminateParticles();
    
    void removeEmptyParticles();
	
	void preAllocateNodes();
	
	void findNodes( const Sorter& buckets, const VectorXs& x, std::vector< Matrix27x3i >& particle_nodes_x, std::vector< Matrix27x3i >& particle_nodes_y, std::vector< Matrix27x3i >& particle_nodes_z );
    
    void findGaussNodes( const Sorter& buckets, const VectorXs& x, std::vector< Matrix27x3i >& particle_nodes_x, std::vector< Matrix27x3i >& particle_nodes_y, std::vector< Matrix27x3i >& particle_nodes_z );
	
	void findNodesPressure( const Sorter& buckets, const VectorXs& x, std::vector< Matrix27x3i >& particle_nodes_p );
	
	void findSolidPhiNodes( const Sorter& buckets, const VectorXs& x, std::vector< Matrix27x3i >& particle_nodes_sphi );
    
    void findEdgeNodes( const Sorter& buckets, const VectorXs& x );
	
	void generateNodes();
	
	void connectPressureNodes();
	
	void connectSolidPhiNodes();
    
    void connectEdgeNodes();
	
	void postAllocateNodes();
	
	void updateParticleDiv();
    
    scalar getMaxVelocity() const;
    
    scalar getMaxFluidVelocity() const;
	
	scalar getDragCoeff( const scalar& psi, const scalar& sat, const scalar& dv, int material ) const;
    
    scalar getPlanarDragCoeff( const scalar& psi, const scalar& sat, const scalar& dv, int material ) const;
	
	scalar getVerticalDiffusivity( const scalar& psi, int material ) const;
    
    scalar getDragCoeffWithOrientation( const scalar& psi, const scalar& sat, const scalar& dv, const Vector3s& orientation, const scalar& shape_factor, int index, int material ) const;
	
	const Vector3s& getBucketMinCorner() const;
	
	scalar getBucketLength() const;
	
	int getNumColors() const;
    
    const std::vector< int >& getParticleToSurfels() const;
	
	const std::vector< VectorXi >& getNodeColorP() const;
	
	const std::vector< VectorXs >& getNodeSurfTensionP() const;
    
    void updateIntersection();
    
    const std::vector< std::vector<RayTriInfo> >& getIntersections() const;
    
    const std::vector<Vector3s>& getFaceWeights() const;
    
    scalar interpolateBucketFluidVelX(const Vector3s& pos) const;
    
    scalar interpolateBucketFluidVelY(const Vector3s& pos) const;
    
    scalar interpolateBucketFluidVelZ(const Vector3s& pos) const;
    
    scalar interpolateBucketSavedFluidVelX(const Vector3s& pos) const;
    
    scalar interpolateBucketSavedFluidVelY(const Vector3s& pos) const;
    
    scalar interpolateBucketSavedFluidVelZ(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidVelX(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidVelY(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidVelZ(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidWeightX(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidWeightY(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidWeightZ(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidPhi(const Vector3s& pos) const;
    
    scalar interpolateBucketLiquidPhi(const Vector3s& pos) const;
    
    scalar interpolateBucketPressure(const Vector3s& pos) const;
    
    scalar interpolateBucketSolidPhiGrad(const Vector3s& pos, Vector3s& grad) const;
    
    scalar interpolateValue(const Vector3s& pos, const std::vector< VectorXs >& phi, const std::vector< VectorXi >& cpidx_phi, const Vector3s& phi_ori, const scalar& default_val);
    
    inline Vector3s nodePosFromBucket(int bucket_idx, int raw_node_idx, const Vector3s& offset) const;
    
    void markInsideOut();
    
    bool isSoft(int pidx) const;
    
    bool isOutsideFluid(int pidx) const;
    
    void expandFluidNodesMarked(int layers);
    
    void updateVelocityDifference();
    
    void saveFluidVelocity();
    
    void saveParticleVelocity();
    
    void updateShapeFactor();
    
    void updateOrientation();
    
    scalar totalFluidVolumeParticles() const;
    
    scalar totalFluidVolumeSoftElasto() const;
	
	void computeDDA();

	void insertSolveGroup(const VectorXi& group);

	const std::vector<VectorXi>& getSolveGroup() const;
	
	void constrainLiquidVelocity();
	
	void updateStrandParamViscosity(const scalar& dt);
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
private:
    int step_count;
    VectorXs m_x; //particle pos
	VectorXs m_rest_x; //particle rest pos
	
    VectorXs m_v; //particle velocity
    VectorXs m_saved_v;
	VectorXs m_dv;
    VectorXs m_fluid_v; //fluid velocity
    VectorXs m_m; //particle mass
	VectorXs m_fluid_m;
	VectorXs m_radius; // particle radius
    VectorXs m_vol; //particle volume
	VectorXs m_rest_vol;
    VectorXs m_shape_factor;
	VectorXs m_fluid_vol;
    VectorXuc m_inside;
	VectorXs m_volume_fraction; // elastic particle fraction
	VectorXs m_rest_volume_fraction;
    VectorXs m_orientation;
	std::vector< VectorXs > m_div;
    std::vector< VectorXs > m_sphere_pattern;
    std::vector< ParticleClassifier > m_classifier;
	
	VectorXs m_particle_rest_length;
    VectorXs m_edge_rest_length;
	VectorXs m_particle_rest_area;
    VectorXs m_face_rest_area;
	MatrixXs m_B; // particle B matrix
    MatrixXs m_fB;
	
	VectorXs m_x_gauss;
	VectorXs m_v_gauss;
	VectorXs m_dv_gauss;
    VectorXs m_fluid_v_gauss; //fluid velocity
	VectorXs m_m_gauss;
	VectorXs m_vol_gauss;
	VectorXs m_rest_vol_gauss;
	VectorXs m_radius_gauss;
	VectorXs m_fluid_m_gauss;
	VectorXs m_fluid_vol_gauss;
	VectorXs m_volume_fraction_gauss;
	VectorXs m_rest_volume_fraction_gauss;
	
	std::vector<unsigned char> m_bucket_marked;
	
    std::vector< std::vector<RayTriInfo> > m_ray_tri_gauss;
	
    MatrixXs m_Fe_gauss; // elastic deformation gradient
	MatrixXs m_d_gauss; // material axis
    MatrixXs m_d_old_gauss; //store last deformation information
    MatrixXs m_D_gauss; //original material axis
    MatrixXs m_D_inv_gauss; //inverse of original material axis
    MatrixXs m_dFe_gauss; //dphidFe
	MatrixXs m_grad_gauss; //3m x 3 (or 2) vectors for computing gradient on Gauss
    MatrixXs m_norm_gauss; // normalized material axis with rigid transformation
	
    std::vector<unsigned char> m_fixed;
	std::vector<bool> m_twist;
    MatrixXi m_edges;
	MatrixXi m_faces; //store the face id
	std::vector<Vector3s> m_face_weights;
	
	std::vector<Vector3i> m_face_inv_mapping; // location of the face itself in particle_to_face array.
	std::vector<Vector2i> m_edge_inv_mapping; // location of the edge itself in particle_to_face array.
	
	std::vector<int> m_surfels;
	std::vector<int> m_fluids;
    
    std::vector<Vector3s> m_surfel_norms;
	
	std::vector< std::vector<std::pair<int, scalar> > > m_particle_to_face;
    
    std::vector< std::vector<int> > m_particle_to_edge;
    
    std::vector< int > m_particle_to_surfel;
	std::vector< int > m_gauss_to_parameters;
    std::vector< int > m_face_to_parameters;
    std::vector< int > m_edge_to_parameters;
	
	std::vector< Matrix27x3i > m_particle_nodes_x;
	std::vector< Matrix27x3i > m_particle_nodes_y;
	std::vector< Matrix27x3i > m_particle_nodes_z;
	std::vector< Matrix27x3i > m_particle_nodes_solid_phi;
	std::vector< Matrix27x3i > m_particle_nodes_p;
	
	std::vector< Matrix27x3i > m_gauss_nodes_x;
	std::vector< Matrix27x3i > m_gauss_nodes_y;
	std::vector< Matrix27x3i > m_gauss_nodes_z;
	
	std::vector< Vector27s > m_particle_weights_p;
	std::vector< Matrix27x4s > m_particle_weights;
	
	std::vector< Matrix27x3s > m_particle_grads_x;
	std::vector< Matrix27x3s > m_gauss_grads_x;
	std::vector< Matrix27x3s > m_particle_grads_y;
	std::vector< Matrix27x3s > m_gauss_grads_y;
	std::vector< Matrix27x3s > m_particle_grads_z;
	std::vector< Matrix27x3s > m_gauss_grads_z;
	std::vector< Matrix27x3s > m_particle_grads_solid_phi;
	
	
	// bucket id -> nodes -> pairs of (particle id, id in particle neighborhoods)
	std::vector< std::vector< std::vector< std::pair< int, int > > > > m_node_particles_x;
	std::vector< std::vector< std::vector< std::pair< int, int > > > > m_node_particles_y;
	std::vector< std::vector< std::vector< std::pair< int, int > > > > m_node_particles_z;
	std::vector< std::vector< std::vector< std::pair< int, int > > > > m_node_particles_p;
	
    std::vector< Matrix27x3s > m_gauss_dFdx_x;
    std::vector< Matrix27x3s > m_gauss_dFdx_y;
    std::vector< Matrix27x3s > m_gauss_dFdx_z;
	
	std::vector< int > m_particle_group;
	
	Vector3s m_bbx_min;
	Vector3s m_bbx_max;
	
	scalar m_bucket_size;
	int m_num_colors;
	int m_num_nodes;
	int m_num_bucket_colors;
	Vector3s m_bucket_mincorner;
	Vector3s m_grid_mincorner;
	
	Sorter m_particle_buckets;
	Sorter m_gauss_buckets;
	Sorter m_particle_cells;
	
	std::vector< VectorXi > m_node_cpidx_x_tmp;
	std::vector< VectorXi > m_node_cpidx_y_tmp;
	std::vector< VectorXi > m_node_cpidx_z_tmp;
	std::vector< VectorXi > m_node_cpidx_p_tmp;
	std::vector< VectorXi > m_node_cpidx_solid_phi_tmp;
	
	std::vector< VectorXi > m_node_cpidx_ex_tmp;
	std::vector< VectorXi > m_node_cpidx_ey_tmp;
	std::vector< VectorXi > m_node_cpidx_ez_tmp;
	
	std::vector< VectorXi > m_node_cpidx_x; // bucket id -> node to compressed index
	std::vector< VectorXi > m_node_cpidx_y; // bucket id -> node to compressed index
	std::vector< VectorXi > m_node_cpidx_z; // bucket id -> node to compressed index
	
	std::vector< VectorXi > m_node_cpidx_p; // bucket id -> node to compressed index of pressure
	
	std::vector< VectorXi > m_node_cpidx_solid_phi; // bucket id -> node to compressed index of solid phi
	
    std::vector< VectorXi > m_node_cpidx_ex; // bucket id -> node to compressed index of edge nodes
    std::vector< VectorXi > m_node_cpidx_ey; // bucket id -> node to compressed index of edge nodes
    std::vector< VectorXi > m_node_cpidx_ez; // bucket id -> node to compressed index of edge nodes
    
	std::vector< VectorXs > m_node_pos_x; // bucket id -> nodes pos
	std::vector< VectorXs > m_node_pos_y; // bucket id -> nodes pos
	std::vector< VectorXs > m_node_pos_z;
	std::vector< VectorXs > m_node_pos_p;
    
    std::vector< VectorXs > m_node_pos_ex; // bucket id -> nodes pos
    std::vector< VectorXs > m_node_pos_ey; // bucket id -> nodes pos
    std::vector< VectorXs > m_node_pos_ez;
    
    std::vector< VectorXi > m_node_indices_x;
    std::vector< VectorXi > m_node_indices_y;
    std::vector< VectorXi > m_node_indices_z;
    std::vector< VectorXi > m_node_indices_p;
	
	std::vector< VectorXs > m_node_pos_solid_phi;
    
    std::vector< VectorXs > m_node_shape_factor_x;
    std::vector< VectorXs > m_node_shape_factor_y;
    std::vector< VectorXs > m_node_shape_factor_z;
	
	std::vector< VectorXi > m_node_pressure_neighbors; // bucket id -> 6x2 neighbors of pressure nodes
	std::vector< VectorXi > m_node_pp_neighbors; // bucket id -> 18x2 pressure neighbors of pressure nodes
	
	std::vector< VectorXi > m_node_index_pressure_x; // bucket id -> 2x2 neighbors (left/right) to the pressure nodes
	std::vector< VectorXi > m_node_index_pressure_y; // bucket id -> 2x2 neighbors (down/up) to the pressure nodes
	std::vector< VectorXi > m_node_index_pressure_z; // bucket id -> 2x2 neighbors (back/front) to the pressure nodes
	
	std::vector< VectorXi > m_node_index_solid_phi_x; // bucket id -> 4x2 neighbors to the solid phi nodes
	std::vector< VectorXi > m_node_index_solid_phi_y; // bucket id -> 4x2 neighbors to the solid phi nodes
	std::vector< VectorXi > m_node_index_solid_phi_z; // bucket id -> 4x2 neighbors to the solid phi nodes
    
    std::vector< VectorXi > m_node_index_edge_x; // bucket id -> 4x2 neighbors to the edge nodes (ey, ez)
    std::vector< VectorXi > m_node_index_edge_y; // bucket id -> 4x2 neighbors to the edge nodes (ex, ez)
    std::vector< VectorXi > m_node_index_edge_z; // bucket id -> 4x2 neighbors to the edge nodes (ex, ey)
	
	std::vector< VectorXs > m_node_solid_phi;
	std::vector< VectorXs > m_node_liquid_phi;
	std::vector< VectorXs > m_node_pressure;
	std::vector< VectorXs > m_node_cell_solid_phi;
	
	std::vector< VectorXuc > m_node_state_u;
	std::vector< VectorXuc > m_node_state_v;
	std::vector< VectorXuc > m_node_state_w;
	
    // volume fractions for viscosity
    std::vector< VectorXs > m_node_liquid_c_vf;
    std::vector< VectorXs > m_node_liquid_u_vf;
    std::vector< VectorXs > m_node_liquid_v_vf;
    std::vector< VectorXs > m_node_liquid_w_vf;
    std::vector< VectorXs > m_node_liquid_ex_vf;
    std::vector< VectorXs > m_node_liquid_ey_vf;
    std::vector< VectorXs > m_node_liquid_ez_vf;
    
    std::vector< VectorXs > m_node_orientation_x;
    std::vector< VectorXs > m_node_orientation_y;
    std::vector< VectorXs > m_node_orientation_z;
	
	std::vector< VectorXs > m_node_sat_x; // bucket id -> nodes saturation
	std::vector< VectorXs > m_node_sat_y; // bucket id -> nodes saturation
	std::vector< VectorXs > m_node_sat_z; // bucket id -> nodes saturation
	
	std::vector< VectorXs > m_node_sat_p; // bucket id -> nodes saturation
	
	std::vector< VectorXs > m_node_pore_pressure_p; // bucket id -> nodes saturation
	
	std::vector< VectorXs > m_node_psi_x;
	std::vector< VectorXs > m_node_psi_y;
	std::vector< VectorXs > m_node_psi_z;
	
	std::vector< VectorXs > m_node_psi_p;

	std::vector< VectorXs > m_node_solid_vel_x;
	std::vector< VectorXs > m_node_solid_vel_y;
	std::vector< VectorXs > m_node_solid_vel_z;
	
	std::vector< VectorXuc > m_node_liquid_valid_x;
	std::vector< VectorXuc > m_node_liquid_valid_y;
	std::vector< VectorXuc > m_node_liquid_valid_z;
	
	std::vector< VectorXs > m_node_solid_weight_x; // bucket id -> node solid weight
	std::vector< VectorXs > m_node_solid_weight_y;
	std::vector< VectorXs > m_node_solid_weight_z;
	
	std::vector< VectorXs > m_node_vel_x; // bucket id -> nodes velocity
	std::vector< VectorXs > m_node_vel_y;
	std::vector< VectorXs > m_node_vel_z;
	
	std::vector< VectorXs > m_node_mass_x; // bucket id -> nodes mass
	std::vector< VectorXs > m_node_mass_y;
	std::vector< VectorXs > m_node_mass_z;
	
	std::vector< VectorXs > m_node_raw_weight_x; // bucket id -> nodes volume
	std::vector< VectorXs > m_node_raw_weight_y;
	std::vector< VectorXs > m_node_raw_weight_z;
	
	std::vector< VectorXs > m_node_vol_fluid_x; // bucket id -> nodes volume
	std::vector< VectorXs > m_node_vol_fluid_y;
	std::vector< VectorXs > m_node_vol_fluid_z;
	
	std::vector< VectorXs > m_node_vol_pure_fluid_x; // bucket id -> nodes volume
	std::vector< VectorXs > m_node_vol_pure_fluid_y;
	std::vector< VectorXs > m_node_vol_pure_fluid_z;
	
	std::vector< VectorXs > m_node_vol_x; // bucket id -> nodes volume
	std::vector< VectorXs > m_node_vol_y;
	std::vector< VectorXs > m_node_vol_z;
	
	std::vector< VectorXs > m_node_vel_fluid_x; // bucket id -> nodes velocity
	std::vector< VectorXs > m_node_vel_fluid_y;
	std::vector< VectorXs > m_node_vel_fluid_z;
    
    std::vector< VectorXs > m_node_vel_saved_fluid_x; // bucket id -> nodes velocity
    std::vector< VectorXs > m_node_vel_saved_fluid_y;
    std::vector< VectorXs > m_node_vel_saved_fluid_z;
	
	std::vector< VectorXs > m_node_mass_fluid_x; // bucket id -> nodes mass
	std::vector< VectorXs > m_node_mass_fluid_y;
	std::vector< VectorXs > m_node_mass_fluid_z;
	
	std::vector< std::shared_ptr<StrandForce> > m_strands;
    std::vector<bool> m_is_strand_tip;

    std::vector< VectorXi > m_solve_groups;
    
    std::vector< MatrixXi > m_gauss_bucket_neighbors;
	
	LiquidInfo m_liquid_info;
	
    // Forces. Note that the scene inherits responsibility for deleting forces.
    std::vector< std::shared_ptr<Force> > m_forces;
    
    std::vector< std::shared_ptr<AttachForce> > m_attach_forces;

	std::vector< std::shared_ptr<StrandParameters> > m_strandParameters;
	
	std::vector< Vector3s > m_group_pos;
	std::vector< Eigen::Quaternion<scalar> > m_group_rot;
	
	std::vector< Vector3s > m_group_prev_pos;
	std::vector< Eigen::Quaternion<scalar> > m_group_prev_rot;
	
	std::vector< std::shared_ptr<Script> > m_scripts;
    
    std::vector< scalar > m_shooting_vol_accum;
	
	std::vector< std::shared_ptr< DistanceField > > m_group_distance_field;
	
	std::vector< std::shared_ptr< DistanceField > > m_distance_fields;
};

#endif
