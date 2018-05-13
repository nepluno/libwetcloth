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

#ifndef __LINEARIZED_IMPLICIT_EULER__
#define __LINEARIZED_IMPLICIT_EULER__

#include <Eigen/Core>
#include <iostream>

#include "SceneStepper.h"
#include "MathUtilities.h"
#include "StringUtilities.h"
#include "array3.h"
#include "pcgsolver/sparse_matrix.h"

class LinearizedImplicitEuler : public SceneStepper
{
public:
	LinearizedImplicitEuler( const scalar& criterion, const scalar& pressure_criterion, int maxiters, int manifold_substeps, int viscosity_substeps );
	
	virtual ~LinearizedImplicitEuler();
	
	virtual bool stepScene( TwoDScene& scene, scalar dt );
	
	virtual bool stepVelocity( TwoDScene& scene, scalar dt );
    
    virtual bool projectFine( TwoDScene& scene, scalar dt );

	virtual bool acceptVelocity( TwoDScene& scene );
    
    virtual bool stepImplicitElasto( TwoDScene& scene, scalar dt );
    
    virtual bool stepImplicitElastoDiagonalPCR( TwoDScene& scene, scalar dt );
    
    virtual bool stepImplicitElastoDiagonalPCG( TwoDScene& scene, scalar dt );
    
    virtual bool stepImplicitElastoAMGPCG( TwoDScene& scene, scalar dt );
    
    virtual bool applyPressureDragElasto( TwoDScene& scene, scalar dt );
    
    virtual bool applyPressureDragFluid( TwoDScene& scene, scalar dt );
    
    virtual bool solveBiCGSTAB( TwoDScene& scene, scalar dt );
	
	virtual bool manifoldPropagate( TwoDScene& scene, scalar dt );
    
    virtual scalar computeDivergence( TwoDScene& scene );
    
    virtual void pushFluidVelocity();
    
    virtual void popFluidVelocity();
    
    virtual void pushElastoVelocity();
    
    virtual void popElastoVelocity();
	
	virtual std::string getName() const;
	
private:
	void zeroFixedDoFs( const TwoDScene& scene, VectorXs& vec );
		
	void performLocalSolve( const TwoDScene& scene,
						   const std::vector< VectorXs >& node_rhs_x,
						   const std::vector< VectorXs >& node_rhs_y,
						   const std::vector< VectorXs >& node_rhs_z,
						   const std::vector< VectorXs >& node_mass_x,
						   const std::vector< VectorXs >& node_mass_y,
						   const std::vector< VectorXs >& node_mass_z,
						   std::vector< VectorXs >& out_node_vec_x,
						   std::vector< VectorXs >& out_node_vec_y,
						   std::vector< VectorXs >& out_node_vec_z );
	
	void performInvLocalSolve( const TwoDScene& scene,
						   const std::vector< VectorXs >& node_rhs_x,
						   const std::vector< VectorXs >& node_rhs_y,
						   const std::vector< VectorXs >& node_rhs_z,
						   const std::vector< VectorXs >& node_inv_mass_x,
						   const std::vector< VectorXs >& node_inv_mass_y,
						   const std::vector< VectorXs >& node_inv_mass_z,
						   std::vector< VectorXs >& out_node_vec_x,
						   std::vector< VectorXs >& out_node_vec_y,
						   std::vector< VectorXs >& out_node_vec_z );
	
	void performGlobalMultiply( const TwoDScene& scene, const scalar& dt,
							   const std::vector< VectorXs >& node_m_x,
							   const std::vector< VectorXs >& node_m_y,
							   const std::vector< VectorXs >& node_m_z,
							   const std::vector< VectorXs >& node_v_x,
							   const std::vector< VectorXs >& node_v_y,
							   const std::vector< VectorXs >& node_v_z,
							   std::vector< VectorXs >& out_node_vec_x,
							   std::vector< VectorXs >& out_node_vec_y,
							   std::vector< VectorXs >& out_node_vec_z );
    
    void performGlobalMultiplyBiCGSTAB( const TwoDScene& scene, const scalar& dt,
                                       const std::vector< VectorXs >& node_v_s_x,
                                       const std::vector< VectorXs >& node_v_s_y,
                                       const std::vector< VectorXs >& node_v_s_z,
                                       const std::vector< VectorXs >& node_v_f_x,
                                       const std::vector< VectorXs >& node_v_f_y,
                                       const std::vector< VectorXs >& node_v_f_z,
                                       const std::vector< VectorXs >& node_v_p,
                                       std::vector< VectorXs >& out_node_vec_s_x,
                                       std::vector< VectorXs >& out_node_vec_s_y,
                                       std::vector< VectorXs >& out_node_vec_s_z,
                                       std::vector< VectorXs >& out_node_vec_f_x,
                                       std::vector< VectorXs >& out_node_vec_f_y,
                                       std::vector< VectorXs >& out_node_vec_f_z,
                                       std::vector< VectorXs >& out_node_vec_p);
	
    void constructNodeForceCoarse( TwoDScene& scene, const scalar& dt,
                                  Array3s& node_rhs_fluid_x,
                                  Array3s& node_rhs_fluid_y,
                                  Array3s& node_rhs_fluid_z );
    
    void performInvLocalSolveCoarse( TwoDScene& scene, 
                                    const Array3s& node_rhs_x,
                                    const Array3s& node_rhs_y,
                                    const Array3s& node_rhs_z,
                                    const Array3s& node_inv_mass_x,
                                    const Array3s& node_inv_mass_y,
                                    const Array3s& node_inv_mass_z,
                                    Array3s& out_node_vec_x,
                                    Array3s& out_node_vec_y,
                                    Array3s& out_node_vec_z );
    
    void constructHDVCoarse( TwoDScene& scene, const scalar& dt,
                            Array3s& node_hdv_x,
                            Array3s& node_hdv_y,
                            Array3s& node_hdv_z,
                            Array3s& node_hdvs_x,
                            Array3s& node_hdvs_y,
                            Array3s& node_hdvs_z,
                            const Array3s& node_v_x,
                            const Array3s& node_v_y,
                            const Array3s& node_v_z,
                            const Array3s& node_fluid_v_x,
                            const Array3s& node_fluid_v_y,
                            const Array3s& node_fluid_v_z);
    
    void constructInvMDVCoarse( TwoDScene& scene,
                         Array3s& node_inv_mdv_x,
                         Array3s& node_inv_mdv_y,
                         Array3s& node_inv_mdv_z,
                         const Array3s& node_hdv_x,
                         const Array3s& node_hdv_y,
                         const Array3s& node_hdv_z);
    
    void constructMsDVsCoarse( TwoDScene& scene,
                        Array3s& node_msdv2_x,
                        Array3s& node_msdv2_y,
                        Array3s& node_msdv2_z,
                        Array3s& node_inv_msdv2_x,
                        Array3s& node_inv_msdv2_y,
                        Array3s& node_inv_msdv2_z,
                        const Array3s& node_hdvs_x,
                        const Array3s& node_hdvs_y,
                        const Array3s& node_hdvs_z);
    
    void constructPsiSFCoarse( TwoDScene& scene,
                              Array3s& node_psi_sf_x,
                              Array3s& node_psi_sf_y,
                              Array3s& node_psi_sf_z,
                              const Array3s& node_inv_mdv_x,
                              const Array3s& node_inv_mdv_y,
                              const Array3s& node_inv_mdv_z,
                              const Array3s& node_hdv_x,
                              const Array3s& node_hdv_y,
                              const Array3s& node_hdv_z );
    
    void constructPsiFSCoarse( TwoDScene& scene,
                              Array3s& node_psi_fs_x,
                              Array3s& node_psi_fs_y,
                              Array3s& node_psi_fs_z,
                              const Array3s& node_inv_mdvs_x,
                              const Array3s& node_inv_mdvs_y,
                              const Array3s& node_inv_mdvs_z,
                              const Array3s& node_hdvs_x,
                              const Array3s& node_hdvs_y,
                              const Array3s& node_hdvs_z );
    
	void constructNodeForce( TwoDScene& scene, const scalar& dt, std::vector< VectorXs >& node_rhs_x, std::vector< VectorXs >& node_rhs_y, std::vector< VectorXs >& node_rhs_z, std::vector< VectorXs >& node_rhs_fluid_x, std::vector< VectorXs >& node_rhs_fluid_y, std::vector< VectorXs >& node_rhs_fluid_z );
	
    void addSolidDragRHS( TwoDScene& scene,
                         const std::vector< VectorXs >& node_vel_x,
                         const std::vector< VectorXs >& node_vel_y,
                         const std::vector< VectorXs >& node_vel_z,
                         std::vector< VectorXs >& node_rhs_x,
                         std::vector< VectorXs >& node_rhs_y,
                         std::vector< VectorXs >& node_rhs_z );
    
	void addFluidDragRHS( TwoDScene& scene,
                         const std::vector< VectorXs >& node_fluid_vel_x,
                         const std::vector< VectorXs >& node_fluid_vel_y,
                         const std::vector< VectorXs >& node_fluid_vel_z,
                         std::vector< VectorXs >& node_rhs_x,
                         std::vector< VectorXs >& node_rhs_y,
                         std::vector< VectorXs >& node_rhs_z );
	
	void addSolidDrag( TwoDScene& scene,
                      const std::vector< VectorXs >& node_vel_x,
                      const std::vector< VectorXs >& node_vel_y,
                      const std::vector< VectorXs >& node_vel_z,
                      std::vector< VectorXs >& node_fluid_vel_x,
                      std::vector< VectorXs >& node_fluid_vel_y,
                      std::vector< VectorXs >& node_fluid_vel_z  );
	
	void constructHDV( TwoDScene& scene, const scalar& dt );
	
	void constructInvMDV( TwoDScene& scene);
	
	void constructMsDVs( TwoDScene& scene);
	
	void constructPsiSF( TwoDScene& scene );
	
	void constructHessianPreProcess( TwoDScene& scene, const scalar& dt );
    
    void constructHessianPostProcess( TwoDScene& scene, const scalar& dt );
	
//    std::vector< Eigen::SimplicialLDLT< SparseXs >* > m_local_solvers;
	
//    SparseXs m_A;
    std::vector< std::pair<int, int> > m_triA_sup;
    TripletXs m_triA;
    VectorXs m_multiply_buffer;
    VectorXs m_pre_mult_buffer;
	
	std::vector< std::vector< std::vector< Matrix3s > > > m_gauss_ddphidfdf;
	
	const scalar m_pcg_criterion;
    const scalar m_pressure_criterion;
	const int m_maxiters;
	const int m_manifold_substeps;
	const int m_viscosity_substeps;
	
	std::vector< VectorXs > m_node_rhs_x;
	std::vector< VectorXs > m_node_rhs_y;
	std::vector< VectorXs > m_node_rhs_z;
	std::vector< VectorXs > m_node_v_plus_x;
	std::vector< VectorXs > m_node_v_plus_y;
	std::vector< VectorXs > m_node_v_plus_z;
	std::vector< VectorXs > m_node_v_fluid_plus_x;
	std::vector< VectorXs > m_node_v_fluid_plus_y;
	std::vector< VectorXs > m_node_v_fluid_plus_z;
    std::vector< VectorXs > m_node_v_0_x;
    std::vector< VectorXs > m_node_v_0_y;
    std::vector< VectorXs > m_node_v_0_z;
    std::vector< VectorXs > m_node_v_tmp_x;
    std::vector< VectorXs > m_node_v_tmp_y;
    std::vector< VectorXs > m_node_v_tmp_z;
	std::vector< VectorXs > m_node_r_x; //r0
	std::vector< VectorXs > m_node_r_y;
	std::vector< VectorXs > m_node_r_z;
	std::vector< VectorXs > m_node_z_x; // s
	std::vector< VectorXs > m_node_z_y;
	std::vector< VectorXs > m_node_z_z;
	std::vector< VectorXs > m_node_p_x; // p
	std::vector< VectorXs > m_node_p_y;
	std::vector< VectorXs > m_node_p_z;
	std::vector< VectorXs > m_node_q_x; // h
	std::vector< VectorXs > m_node_q_y;
	std::vector< VectorXs > m_node_q_z;
    std::vector< VectorXs > m_node_w_x; // v
    std::vector< VectorXs > m_node_w_y;
    std::vector< VectorXs > m_node_w_z;
    std::vector< VectorXs > m_node_t_x; // t
    std::vector< VectorXs > m_node_t_y;
    std::vector< VectorXs > m_node_t_z;
	
	std::vector< VectorXs > m_node_rhs_p;
	std::vector< VectorXs > m_node_r_p;
	std::vector< VectorXs > m_node_z_p;
	std::vector< VectorXs > m_node_p_p;
	std::vector< VectorXs > m_node_q_p;
	std::vector< VectorXs > m_node_jacobi_precond;
	
	std::vector< VectorXs > m_node_rhs_sat;
	std::vector< VectorXs > m_node_perm_sat_x;
	std::vector< VectorXs > m_node_perm_sat_y;
	std::vector< VectorXs > m_node_perm_sat_z;
	std::vector< VectorXs > m_node_r_sat;
	std::vector< VectorXs > m_node_z_sat;
	std::vector< VectorXs > m_node_p_sat;
	std::vector< VectorXs > m_node_q_sat;
	std::vector< VectorXs > m_node_jacobi_precond_sat;
	
	std::vector< VectorXs > m_node_rhs_fluid_x;
	std::vector< VectorXs > m_node_rhs_fluid_y;
	std::vector< VectorXs > m_node_rhs_fluid_z;
    
    std::vector< VectorXs > m_node_hdvm_x; // hDVm
    std::vector< VectorXs > m_node_hdvm_y;
    std::vector< VectorXs > m_node_hdvm_z;
    
    std::vector< VectorXs > m_node_epsilon_x; // 1 - psi
    std::vector< VectorXs > m_node_epsilon_y;
    std::vector< VectorXs > m_node_epsilon_z;
    
    std::vector< VectorXs > m_node_mshdvm_x; // M_s + hDVm
    std::vector< VectorXs > m_node_mshdvm_y;
    std::vector< VectorXs > m_node_mshdvm_z;
    
    std::vector< VectorXs > m_node_damped_x; // damped M_s
    std::vector< VectorXs > m_node_damped_y;
    std::vector< VectorXs > m_node_damped_z;
    
    std::vector< VectorXs > m_node_mfhdvm_x; // M_f + hDVm
    std::vector< VectorXs > m_node_mfhdvm_y;
    std::vector< VectorXs > m_node_mfhdvm_z;
    
    std::vector< VectorXs > m_node_inv_mfhdvm_x; // (M_f+hDVm)^{-1}
    std::vector< VectorXs > m_node_inv_mfhdvm_y;
    std::vector< VectorXs > m_node_inv_mfhdvm_z;
	
    std::vector< VectorXs > m_node_mfhdvm_hdvm_x; // (M_f+hDVm)^{-1}hDVm
	std::vector< VectorXs > m_node_mfhdvm_hdvm_y;
	std::vector< VectorXs > m_node_mfhdvm_hdvm_z;
	
	std::vector< VectorXs > m_node_mshdvm_hdvm_x; // (M_s+hDVm)^{-1}hDVm
	std::vector< VectorXs > m_node_mshdvm_hdvm_y;
	std::vector< VectorXs > m_node_mshdvm_hdvm_z;
	
    std::vector< VectorXs > m_node_inv_C_x; // (M_f+[hdv]*M_s)^{-1}
	std::vector< VectorXs > m_node_inv_C_y;
	std::vector< VectorXs > m_node_inv_C_z;
	
	std::vector< VectorXs > m_node_Cs_x; // M_s+[hdvs]*M_f
	std::vector< VectorXs > m_node_Cs_y;
	std::vector< VectorXs > m_node_Cs_z;
	
    std::vector< VectorXs > m_node_inv_Cs_x; // (M_s+[hdvs]*M_f)^{-1}
	std::vector< VectorXs > m_node_inv_Cs_y;
	std::vector< VectorXs > m_node_inv_Cs_z;
	
	std::vector< VectorXs > m_node_psi_sf_x;
	std::vector< VectorXs > m_node_psi_sf_y;
	std::vector< VectorXs > m_node_psi_sf_z;
	
	std::vector< VectorXs > m_node_psi_fs_x;
	std::vector< VectorXs > m_node_psi_fs_y;
	std::vector< VectorXs > m_node_psi_fs_z;
    
    // bicgstab
    std::vector< VectorXs > m_node_bi_r_S_x; //r
    std::vector< VectorXs > m_node_bi_r_S_y;
    std::vector< VectorXs > m_node_bi_r_S_z;
    std::vector< VectorXs > m_node_bi_r_hat_S_x; //r0
    std::vector< VectorXs > m_node_bi_r_hat_S_y;
    std::vector< VectorXs > m_node_bi_r_hat_S_z;
    std::vector< VectorXs > m_node_bi_s_S_x; // s
    std::vector< VectorXs > m_node_bi_s_S_y;
    std::vector< VectorXs > m_node_bi_s_S_z;
    std::vector< VectorXs > m_node_bi_p_S_x; // p
    std::vector< VectorXs > m_node_bi_p_S_y;
    std::vector< VectorXs > m_node_bi_p_S_z;
    std::vector< VectorXs > m_node_bi_h_S_x; // h
    std::vector< VectorXs > m_node_bi_h_S_y;
    std::vector< VectorXs > m_node_bi_h_S_z;
    std::vector< VectorXs > m_node_bi_t_S_x; // h
    std::vector< VectorXs > m_node_bi_t_S_y;
    std::vector< VectorXs > m_node_bi_t_S_z;
    std::vector< VectorXs > m_node_bi_v_S_x; // h
    std::vector< VectorXs > m_node_bi_v_S_y;
    std::vector< VectorXs > m_node_bi_v_S_z;
    
    std::vector< VectorXs > m_node_bi_r_L_x; //r
    std::vector< VectorXs > m_node_bi_r_L_y;
    std::vector< VectorXs > m_node_bi_r_L_z;
    std::vector< VectorXs > m_node_bi_r_hat_L_x; //r0
    std::vector< VectorXs > m_node_bi_r_hat_L_y;
    std::vector< VectorXs > m_node_bi_r_hat_L_z;
    std::vector< VectorXs > m_node_bi_s_L_x; // s
    std::vector< VectorXs > m_node_bi_s_L_y;
    std::vector< VectorXs > m_node_bi_s_L_z;
    std::vector< VectorXs > m_node_bi_p_L_x; // p
    std::vector< VectorXs > m_node_bi_p_L_y;
    std::vector< VectorXs > m_node_bi_p_L_z;
    std::vector< VectorXs > m_node_bi_h_L_x; // h
    std::vector< VectorXs > m_node_bi_h_L_y;
    std::vector< VectorXs > m_node_bi_h_L_z;
    std::vector< VectorXs > m_node_bi_t_L_x; // h
    std::vector< VectorXs > m_node_bi_t_L_y;
    std::vector< VectorXs > m_node_bi_t_L_z;
    std::vector< VectorXs > m_node_bi_v_L_x; // h
    std::vector< VectorXs > m_node_bi_v_L_y;
    std::vector< VectorXs > m_node_bi_v_L_z;
    
    std::vector< VectorXs > m_node_bi_r_P; //r
    std::vector< VectorXs > m_node_bi_r_hat_P; //r0
    std::vector< VectorXs > m_node_bi_s_P; // s
    std::vector< VectorXs > m_node_bi_p_P; // p
    std::vector< VectorXs > m_node_bi_h_P; // h
    std::vector< VectorXs > m_node_bi_t_P; // h
    std::vector< VectorXs > m_node_bi_v_P; // h
    
    std::stack< std::vector< VectorXs > > m_fluid_vel_stack;
    std::stack< std::vector< VectorXs > > m_elasto_vel_stack;
    
    std::vector<double> m_arr_pressure_rhs;
    robertbridson::SparseMatrix<scalar> m_arr_pressure_matrix;

    std::vector<double> m_fine_pressure_rhs;
    robertbridson::SparseMatrix<scalar> m_fine_pressure_matrix;
    std::vector< VectorXi > m_fine_global_indices;
    
    SparseXs m_A;
    std::vector< VectorXi > m_node_global_indices_x;
    std::vector< VectorXi > m_node_global_indices_y;
    std::vector< VectorXi > m_node_global_indices_z;
    std::vector< Vector3i > m_effective_node_indices;
    std::vector< Vector3i > m_dof_ijk;
    TripletXs m_tri_W;
    SparseXs m_W;
    TripletXs m_tri_M;
    SparseXs m_M;
    std::vector< double > m_elasto_rhs;
    std::vector< double > m_elasto_result;
    robertbridson::SparseMatrix<scalar> m_H;
};

#endif
