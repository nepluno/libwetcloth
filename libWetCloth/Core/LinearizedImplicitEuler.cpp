//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "LinearizedImplicitEuler.h"
#include "ThreadUtils.h"
#include "MathUtilities.h"
#include "sorter.h"
#include "Pressure.h"
#include "Viscosity.h"
#include "array3_utils.h"
#include "AlgebraicMultigrid.h"

#include <unordered_map>

//#define OPTIMIZE_SAT
//#define CHECK_EQU_24

LinearizedImplicitEuler::LinearizedImplicitEuler(const scalar& criterion, const scalar& pressure_criterion, const scalar& quasi_static_criterion, const scalar& viscous_criterion, int maxiters, int manifold_substeps, int viscosity_substeps, int surf_tension_substeps)
	: SceneStepper(), m_pcg_criterion(criterion), m_pressure_criterion(pressure_criterion), m_quasi_static_criterion(quasi_static_criterion), m_viscous_criterion(viscous_criterion), m_maxiters(maxiters), m_manifold_substeps(manifold_substeps), m_viscosity_substeps(viscosity_substeps), m_surf_tension_substeps(surf_tension_substeps)
{}

LinearizedImplicitEuler::~LinearizedImplicitEuler()
{
}

bool LinearizedImplicitEuler::advectSurfTension( TwoDScene& scene, scalar dt )
{
	if (!scene.useSurfTension()) return true;

	const scalar subdt = dt / (scalar) m_surf_tension_substeps;

	for (int i = 0; i < m_surf_tension_substeps; ++i)
	{
		scene.updateColorP();
		scene.renormalizeLiquidPhi();
		scene.updateCurvatureP();
		scene.advectCurvatureP(subdt);

		std::cout << "[surface tension advect: " << (i + 1) << " / " << m_surf_tension_substeps << "]" << std::endl;
	}

	return true;
}

bool LinearizedImplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
	return true;
}

void LinearizedImplicitEuler::performInvLocalSolve( const TwoDScene& scene,
        const std::vector< VectorXs >& node_rhs_x,
        const std::vector< VectorXs >& node_rhs_y,
        const std::vector< VectorXs >& node_rhs_z,
        const std::vector< VectorXs >& node_inv_mass_x,
        const std::vector< VectorXs >& node_inv_mass_y,
        const std::vector< VectorXs >& node_inv_mass_z,
        std::vector< VectorXs >& out_node_vec_x,
        std::vector< VectorXs >& out_node_vec_y,
        std::vector< VectorXs >& out_node_vec_z )
{
	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		const VectorXs& bucket_node_masses_x = node_inv_mass_x[bucket_idx];
		const VectorXs& bucket_node_masses_y = node_inv_mass_y[bucket_idx];
		const VectorXs& bucket_node_masses_z = node_inv_mass_z[bucket_idx];

		const VectorXs& bucket_node_rhs_x = node_rhs_x[bucket_idx];
		const VectorXs& bucket_node_rhs_y = node_rhs_y[bucket_idx];
		const VectorXs& bucket_node_rhs_z = node_rhs_z[bucket_idx];

		VectorXs& bucket_out_node_vec_x = out_node_vec_x[bucket_idx];
		VectorXs& bucket_out_node_vec_y = out_node_vec_y[bucket_idx];
		VectorXs& bucket_out_node_vec_z = out_node_vec_z[bucket_idx];

		for (int i = 0; i < num_nodes; ++i)
		{
			bucket_out_node_vec_x[i] = bucket_node_rhs_x[i] * bucket_node_masses_x[i];
			bucket_out_node_vec_y[i] = bucket_node_rhs_y[i] * bucket_node_masses_y[i];
			bucket_out_node_vec_z[i] = bucket_node_rhs_z[i] * bucket_node_masses_z[i];
		}

		assert(!std::isnan(bucket_out_node_vec_x.sum()));
		assert(!std::isnan(bucket_out_node_vec_y.sum()));
		assert(!std::isnan(bucket_out_node_vec_z.sum()));
	});
}

void LinearizedImplicitEuler::performLocalSolve( const TwoDScene& scene, const VectorXs& rhs, const VectorXs& m, VectorXs& out)
{
	const int num_elasto = scene.getNumSoftElastoParticles();

	threadutils::for_each(0, num_elasto * 4, [&](int pidx) {
		out(pidx) = rhs(pidx) / std::max(1e-20, m(pidx));
	});
}

void LinearizedImplicitEuler::performLocalSolveTwist( const TwoDScene& scene, const VectorXs& rhs, const VectorXs& m, VectorXs& out)
{
	const int num_elasto = scene.getNumSoftElastoParticles();

	threadutils::for_each(0, num_elasto, [&](int pidx) {
		out(pidx) = rhs(pidx) / std::max(1e-20, m(pidx * 4 + 3));
	});
}

void LinearizedImplicitEuler::performLocalSolve( const TwoDScene& scene,
        const std::vector< VectorXs >& node_rhs_x,
        const std::vector< VectorXs >& node_rhs_y,
        const std::vector< VectorXs >& node_rhs_z,
        const std::vector< VectorXs >& node_masses_x,
        const std::vector< VectorXs >& node_masses_y,
        const std::vector< VectorXs >& node_masses_z,
        std::vector< VectorXs >& out_node_vec_x,
        std::vector< VectorXs >& out_node_vec_y,
        std::vector< VectorXs >& out_node_vec_z )
{
	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		const VectorXs& bucket_node_masses_x = node_masses_x[bucket_idx];
		const VectorXs& bucket_node_masses_y = node_masses_y[bucket_idx];
		const VectorXs& bucket_node_masses_z = node_masses_z[bucket_idx];

		const VectorXs& bucket_node_rhs_x = node_rhs_x[bucket_idx];
		const VectorXs& bucket_node_rhs_y = node_rhs_y[bucket_idx];
		const VectorXs& bucket_node_rhs_z = node_rhs_z[bucket_idx];

		VectorXs& bucket_out_node_vec_x = out_node_vec_x[bucket_idx];
		VectorXs& bucket_out_node_vec_y = out_node_vec_y[bucket_idx];
		VectorXs& bucket_out_node_vec_z = out_node_vec_z[bucket_idx];

		for (int i = 0; i < num_nodes; ++i)
		{
			if (bucket_node_masses_x[i] > 1e-20) bucket_out_node_vec_x[i] = bucket_node_rhs_x[i] / bucket_node_masses_x[i];
			else bucket_out_node_vec_x[i] = bucket_node_rhs_x[i];
			if (bucket_node_masses_y[i] > 1e-20) bucket_out_node_vec_y[i] = bucket_node_rhs_y[i] / bucket_node_masses_y[i];
			else bucket_out_node_vec_y[i] = bucket_node_rhs_y[i];
			if (bucket_node_masses_z[i] > 1e-20) bucket_out_node_vec_z[i] = bucket_node_rhs_z[i] / bucket_node_masses_z[i];
			else bucket_out_node_vec_z[i] = bucket_node_rhs_z[i];
		}
	});
}

bool LinearizedImplicitEuler::stepVelocityLagrangian( TwoDScene& scene, scalar dt )
{
	scene.precompute();

	const int num_elasto = scene.getNumSoftElastoParticles();
	const VectorXs& m = scene.getM();
	const VectorXs& v = scene.getV();

	int ndof = num_elasto * 4;

	m_lagrangian_rhs.resize(ndof);
	m_lagrangian_rhs.setZero();

	// accumulate gauss -f
	MatrixXs rhs_gauss(scene.getNumGausses() * 3, 3);
	rhs_gauss.setZero();
	scene.accumulateGaussGradU(rhs_gauss); //force for type 3.

	assert(!std::isnan(rhs_gauss.sum()));

	allocateNodeVectors(scene, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z );
	// map gauss to node
	mapGaussToNode(scene, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z, rhs_gauss);

	// map node to vertices
	mapNodeToSoftParticles(scene, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z, m_lagrangian_rhs);

	// accumulate regular -f
	scene.accumulateGradU(m_lagrangian_rhs);

	assert(!std::isnan(m_lagrangian_rhs.sum()));

	// h * f
	m_lagrangian_rhs *= -dt;

	// h * f + m * v
	m_lagrangian_rhs += VectorXs(m.segment(0, ndof).array() * v.segment(0, ndof).array());

	m_v_plus.resize(ndof);

	performLocalSolve(scene, m_lagrangian_rhs, m, m_v_plus);

	return true;
}

bool LinearizedImplicitEuler::stepVelocity( TwoDScene& scene, scalar dt )
{
	// build node particle pairs
	scene.precompute();

	int ndof_elasto = scene.getNumSoftElastoParticles() * 4;
	const Sorter& buckets = scene.getParticleBuckets();

	allocateNodeVectors(scene, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z );
	allocateNodeVectors(scene, m_node_rhs_fluid_x, m_node_rhs_fluid_y, m_node_rhs_fluid_z );

	allocateNodeVectors(scene, m_node_hdvm_x, m_node_hdvm_y, m_node_hdvm_z );
	allocateNodeVectors(scene, m_node_inv_mfhdvm_x, m_node_inv_mfhdvm_y, m_node_inv_mfhdvm_z );
	allocateNodeVectors(scene, m_node_mfhdvm_hdvm_x, m_node_mfhdvm_hdvm_y, m_node_mfhdvm_hdvm_z );
	allocateNodeVectors(scene, m_node_mshdvm_hdvm_x, m_node_mshdvm_hdvm_y, m_node_mshdvm_hdvm_z );
	allocateNodeVectors(scene, m_node_inv_C_x, m_node_inv_C_y, m_node_inv_C_z );
	allocateNodeVectors(scene, m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z );
	allocateNodeVectors(scene, m_node_Cs_x, m_node_Cs_y, m_node_Cs_z );

	allocateNodeVectors(scene, m_node_psi_sf_x, m_node_psi_sf_y, m_node_psi_sf_z );
	allocateNodeVectors(scene, m_node_psi_fs_x, m_node_psi_fs_y, m_node_psi_fs_z );

	allocateNodeVectors(scene, m_node_damped_x, m_node_damped_y, m_node_damped_z );

	constructNodeForce(scene, dt, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z, m_node_rhs_fluid_x, m_node_rhs_fluid_y, m_node_rhs_fluid_z);

	constructHDV(scene, dt);

	constructInvMDV(scene);

	allocateNodeVectors(scene, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z);

	if (scene.getNumFluidParticles() > 0) {
		const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
		const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
		const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

		performLocalSolve(scene, m_node_rhs_fluid_x, m_node_rhs_fluid_y, m_node_rhs_fluid_z, node_mass_fluid_x, node_mass_fluid_y, node_mass_fluid_z, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z);

		// apply viscosity
		if (scene.getLiquidInfo().compute_viscosity)
		{
			if (scene.getLiquidInfo().implicit_viscosity) {
				stepImplicitViscosityDiagonalPCG(scene, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z,
				                                 m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z, dt);
			} else {
				const scalar sub_dt = dt / (scalar) m_viscosity_substeps;

				for (int i = 0; i < m_viscosity_substeps; ++i)
				{
					m_node_v_tmp_x = m_node_v_fluid_plus_x;
					m_node_v_tmp_y = m_node_v_fluid_plus_y;
					m_node_v_tmp_z = m_node_v_fluid_plus_z;

					viscosity::applyNodeViscosityExplicit(scene, m_node_v_tmp_x, m_node_v_tmp_y, m_node_v_tmp_z, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z, sub_dt);
				}
			}

			// we need recover RHS fluid from new velocity
			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_rhs_fluid_x[bucket_idx] = VectorXs(m_node_v_fluid_plus_x[bucket_idx].array() * node_mass_fluid_x[bucket_idx].array());
				m_node_rhs_fluid_y[bucket_idx] = VectorXs(m_node_v_fluid_plus_y[bucket_idx].array() * node_mass_fluid_y[bucket_idx].array());
				m_node_rhs_fluid_z[bucket_idx] = VectorXs(m_node_v_fluid_plus_z[bucket_idx].array() * node_mass_fluid_z[bucket_idx].array());
			});
		}
	}

	allocateNodeVectors(scene, m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z);

	// construct u_s^*
	if (ndof_elasto > 0 && scene.getLiquidInfo().solve_solid) {
		constructMsDVs(scene);

		const std::vector< VectorXs >& node_mass_x = scene.getNodeMassX();
		const std::vector< VectorXs >& node_mass_y = scene.getNodeMassY();
		const std::vector< VectorXs >& node_mass_z = scene.getNodeMassZ();

		performLocalSolve(scene, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z,
		                  node_mass_x, node_mass_y, node_mass_z,
		                  m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z);

		m_angular_v_plus_buffer.resize(scene.getNumSoftElastoParticles());
		performLocalSolveTwist(scene, m_angular_moment_buffer, scene.getM(), m_angular_v_plus_buffer);
	}

	// record as u_s^* and u_f^*
	acceptVelocity(scene);

	// construct psi_sf and psi_fs
	if (scene.getNumFluidParticles() > 0) {
		constructPsiSF(scene);
	}

	return true;
}


void LinearizedImplicitEuler::performAngularGlobalMultiply( const TwoDScene& scene, const scalar& dt,
        const VectorXs& m,
        const VectorXs& v,
        VectorXs& out)
{
	const int num_elasto = scene.getNumSoftElastoParticles();
	struct Triplets_override
	{
		int m_row, m_col;
		scalar m_value;
	};

	// Ax
	threadutils::for_each(0, num_elasto, [&] (int i) {
		const int idata_start = m_angular_triA_sup[i].first;
		const int idata_end = m_angular_triA_sup[i].second;

		scalar val = 0.0;
		for (int j = idata_start; j < idata_end; ++j)
		{
			const Triplets_override& tri = *((const Triplets_override*) &m_angular_triA[j]);

			val += tri.m_value * v[tri.m_col];
		}
		out[i] = val * dt * dt + m[i * 4 + 3] * v[i];
	});
}

void LinearizedImplicitEuler::performGlobalMultiply( const TwoDScene& scene, const scalar& dt,
        const VectorXs& m,
        const VectorXs& vec,
        VectorXs& out)
{
	const int num_elasto = scene.getNumSoftElastoParticles();

	if (num_elasto == 0) return;

	struct Triplets_override
	{
		int m_row, m_col;
		scalar m_value;
	};

	const int ndof = num_elasto * 4;

	if (m_multiply_buffer.size() != num_elasto * 4) m_multiply_buffer.resize( ndof );
	// Ax
	threadutils::for_each(0, ndof, [&] (int i) {
		const int idata_start = m_triA_sup[i].first;
		const int idata_end = m_triA_sup[i].second;

		scalar val = 0.0;
		for (int j = idata_start; j < idata_end; ++j)
		{
			const Triplets_override& tri = *((const Triplets_override*) &m_triA[j]);

			val += tri.m_value * vec[tri.m_col];
		}
		m_multiply_buffer[i] = val;
	});

	out = m_multiply_buffer * (dt * dt) + VectorXs(m.segment(0, ndof).array() * vec.array());
}

void LinearizedImplicitEuler::performGlobalMultiply( const TwoDScene& scene, const scalar& dt,
        const std::vector< VectorXs >& node_m_x,
        const std::vector< VectorXs >& node_m_y,
        const std::vector< VectorXs >& node_m_z,
        const std::vector< VectorXs >& node_v_x,
        const std::vector< VectorXs >& node_v_y,
        const std::vector< VectorXs >& node_v_z,
        std::vector< VectorXs >& out_node_vec_x,
        std::vector< VectorXs >& out_node_vec_y,
        std::vector< VectorXs >& out_node_vec_z )
{
	const int num_elasto = scene.getNumSoftElastoParticles();

	if (num_elasto == 0) return;
	if (m_multiply_buffer.size() != num_elasto * 4) m_multiply_buffer.resize( num_elasto * 4 );
	if (m_pre_mult_buffer.size() != num_elasto * 4) m_pre_mult_buffer.resize( num_elasto * 4 );

	// Wx
	mapNodeToSoftParticles( scene, node_v_x, node_v_y, node_v_z, m_pre_mult_buffer );

	struct Triplets_override
	{
		int m_row, m_col;
		scalar m_value;
	};

	// AWx
	threadutils::for_each(0, num_elasto * 4, [&] (int i) {
		const int idata_start = m_triA_sup[i].first;
		const int idata_end = m_triA_sup[i].second;

		scalar val = 0.0;
		for (int j = idata_start; j < idata_end; ++j)
		{
			const Triplets_override& tri = *((const Triplets_override*) &m_triA[j]);

			val += tri.m_value * m_pre_mult_buffer[tri.m_col];
		}
		m_multiply_buffer[i] = val;
	});
//    m_multiply_buffer = m_A * m_pre_mult_buffer;

	// W^TAWx
	mapSoftParticlesToNode( scene, out_node_vec_x, out_node_vec_y, out_node_vec_z, m_multiply_buffer );

	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		VectorXs& bucket_node_vec_x = out_node_vec_x[bucket_idx];
		VectorXs& bucket_node_vec_y = out_node_vec_y[bucket_idx];
		VectorXs& bucket_node_vec_z = out_node_vec_z[bucket_idx];

		// (M+h^2(W^TAW+H)x
		bucket_node_vec_x = bucket_node_vec_x * (dt * dt) + VectorXs(node_m_x[bucket_idx].array() * node_v_x[bucket_idx].array());
		bucket_node_vec_y = bucket_node_vec_y * (dt * dt) + VectorXs(node_m_y[bucket_idx].array() * node_v_y[bucket_idx].array());
		bucket_node_vec_z = bucket_node_vec_z * (dt * dt) + VectorXs(node_m_z[bucket_idx].array() * node_v_z[bucket_idx].array());
//		std::cout << bucket_node_vec << std::endl;
	});
}

void LinearizedImplicitEuler::performGlobalMultiply( const TwoDScene& scene, const scalar& dt,
        const std::vector< VectorXs >& node_m_x,
        const std::vector< VectorXs >& node_m_y,
        const std::vector< VectorXs >& node_m_z,
        const std::vector< VectorXs >& node_v_x,
        const std::vector< VectorXs >& node_v_y,
        const std::vector< VectorXs >& node_v_z,
        std::vector< VectorXs >& out_node_vec_x,
        std::vector< VectorXs >& out_node_vec_y,
        std::vector< VectorXs >& out_node_vec_z,
        const VectorXs& m,
        const VectorXs& angular_vec,
        VectorXs& out)
{
	const int num_elasto = scene.getNumSoftElastoParticles();

	if (num_elasto == 0) return;
	if (m_multiply_buffer.size() != num_elasto * 4) m_multiply_buffer.resize( num_elasto * 4 );
	if (m_pre_mult_buffer.size() != num_elasto * 4) m_pre_mult_buffer.resize( num_elasto * 4 );

	// Wx
	mapNodeToSoftParticles( scene, node_v_x, node_v_y, node_v_z, m_pre_mult_buffer );

	// Redistribute angular DOFs
	threadutils::for_each(0, num_elasto, [&] (int i) {
		m_pre_mult_buffer[i * 4 + 3] = angular_vec(i);
	});

	struct Triplets_override
	{
		int m_row, m_col;
		scalar m_value;
	};

	// AWx
	threadutils::for_each(0, num_elasto * 4, [&] (int i) {
		const int idata_start = m_triA_sup[i].first;
		const int idata_end = m_triA_sup[i].second;

		scalar val = 0.0;
		for (int j = idata_start; j < idata_end; ++j)
		{
			const Triplets_override& tri = *((const Triplets_override*) &m_triA[j]);

			val += tri.m_value * m_pre_mult_buffer[tri.m_col];
		}
		m_multiply_buffer[i] = val;
	});


	// grab angular DOFs back
	threadutils::for_each(0, num_elasto, [&] (int i) {
		out(i) = m_multiply_buffer[i * 4 + 3] * (dt * dt) + m(i * 4 + 3) * angular_vec(i);
		m_multiply_buffer[i * 4 + 3] = 0.0;
	});

	// W^TAWx
	mapSoftParticlesToNode( scene, out_node_vec_x, out_node_vec_y, out_node_vec_z, m_multiply_buffer );

	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		VectorXs& bucket_node_vec_x = out_node_vec_x[bucket_idx];
		VectorXs& bucket_node_vec_y = out_node_vec_y[bucket_idx];
		VectorXs& bucket_node_vec_z = out_node_vec_z[bucket_idx];

		// (M+h^2(W^TAW+H)x
		bucket_node_vec_x = bucket_node_vec_x * (dt * dt) + VectorXs(node_m_x[bucket_idx].array() * node_v_x[bucket_idx].array());
		bucket_node_vec_y = bucket_node_vec_y * (dt * dt) + VectorXs(node_m_y[bucket_idx].array() * node_v_y[bucket_idx].array());
		bucket_node_vec_z = bucket_node_vec_z * (dt * dt) + VectorXs(node_m_z[bucket_idx].array() * node_v_z[bucket_idx].array());
//		std::cout << bucket_node_vec << std::endl;
	});
}

void LinearizedImplicitEuler::constructNodeForce( TwoDScene& scene, const scalar& dt, std::vector< VectorXs >& node_rhs_x, std::vector< VectorXs >& node_rhs_y, std::vector< VectorXs >& node_rhs_z, std::vector< VectorXs >& node_rhs_fluid_x, std::vector< VectorXs >& node_rhs_fluid_y, std::vector< VectorXs >& node_rhs_fluid_z )
{
	const int num_elasto = scene.getNumSoftElastoParticles();
	const VectorXs& m = scene.getM();
	const VectorXs& v = scene.getV();

	int ndof = num_elasto * 4;

	VectorXs rhs = VectorXs::Zero(ndof);

	m_angular_moment_buffer.resize(num_elasto);

	if (scene.getLiquidInfo().solve_solid) {
		scene.accumulateGradU(rhs);
		rhs *= -dt;

		assert(!std::isnan(rhs.sum()));

		threadutils::for_each(0, num_elasto, [&] (int pidx) {
			m_angular_moment_buffer[pidx] = rhs[pidx * 4 + 3] + m[pidx * 4 + 3] * v[pidx * 4 + 3];
		});

		mapSoftParticlesToNode(scene, node_rhs_x, node_rhs_y, node_rhs_z, rhs);
	}

	scene.accumulateFluidNodeGradU(node_rhs_fluid_x, node_rhs_fluid_y, node_rhs_fluid_z, -dt);

	////////    // add grad pp to rhs
	pressure::computePorePressureGrads(scene, node_rhs_fluid_x, node_rhs_fluid_y, node_rhs_fluid_z,
	                                   scene.getNodeFluidVolX(), scene.getNodeFluidVolY(), scene.getNodeFluidVolZ(),
	                                   dt);

	if (scene.getLiquidInfo().apply_pore_pressure_solid) {
		pressure::computePorePressureGrads(scene, node_rhs_x, node_rhs_y, node_rhs_z,
		                                   scene.getNodeVolX(), scene.getNodeVolY(), scene.getNodeVolZ(),
		                                   dt);
	}

	if (scene.useSurfTension())
	{
		pressure::applySurfTensionFluid(scene, scene.getNodeSurfTensionP(), node_rhs_fluid_x, node_rhs_fluid_y, node_rhs_fluid_z, scene.getNodeFluidVolX(), scene.getNodeFluidVolY(), scene.getNodeFluidVolZ());
		pressure::applySurfTensionFluid(scene, scene.getNodeSurfTensionP(), node_rhs_x, node_rhs_y, node_rhs_z, scene.getNodeVolX(), scene.getNodeVolY(), scene.getNodeVolZ());
	}

	if (scene.getLiquidInfo().solve_solid) {
		MatrixXs rhs_gauss(scene.getNumGausses() * 3, 3);
		rhs_gauss.setZero();
		scene.accumulateGaussGradU(rhs_gauss); //force for type 3.
		rhs_gauss *= -dt;


		assert(!std::isnan(rhs_gauss.sum()));

		mapGaussToNode(scene, node_rhs_x, node_rhs_y, node_rhs_z, rhs_gauss);
	}
	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (scene.getLiquidInfo().solve_solid) {
			const VectorXs& node_masses_x = scene.getNodeMassX()[bucket_idx];
			const VectorXs& node_masses_y = scene.getNodeMassY()[bucket_idx];
			const VectorXs& node_masses_z = scene.getNodeMassZ()[bucket_idx];

			const VectorXs& node_vel_x = scene.getNodeVelocityX()[bucket_idx];
			const VectorXs& node_vel_y = scene.getNodeVelocityY()[bucket_idx];
			const VectorXs& node_vel_z = scene.getNodeVelocityZ()[bucket_idx];

			node_rhs_x[bucket_idx] += VectorXs(node_masses_x.array() * node_vel_x.array());
			node_rhs_y[bucket_idx] += VectorXs(node_masses_y.array() * node_vel_y.array());
			node_rhs_z[bucket_idx] += VectorXs(node_masses_z.array() * node_vel_z.array());

			assert(!std::isnan(node_rhs_x[bucket_idx].sum()));
			assert(!std::isnan(node_rhs_y[bucket_idx].sum()));
			assert(!std::isnan(node_rhs_z[bucket_idx].sum()));
		}

		const VectorXs& node_masses_fluid_x = scene.getNodeFluidMassX()[bucket_idx];
		const VectorXs& node_masses_fluid_y = scene.getNodeFluidMassY()[bucket_idx];
		const VectorXs& node_masses_fluid_z = scene.getNodeFluidMassZ()[bucket_idx];

		const VectorXs& node_vel_fluid_x = scene.getNodeFluidVelocityX()[bucket_idx];
		const VectorXs& node_vel_fluid_y = scene.getNodeFluidVelocityY()[bucket_idx];
		const VectorXs& node_vel_fluid_z = scene.getNodeFluidVelocityZ()[bucket_idx];

		node_rhs_fluid_x[bucket_idx] += VectorXs(node_masses_fluid_x.array() * node_vel_fluid_x.array());
		node_rhs_fluid_y[bucket_idx] += VectorXs(node_masses_fluid_y.array() * node_vel_fluid_y.array());
		node_rhs_fluid_z[bucket_idx] += VectorXs(node_masses_fluid_z.array() * node_vel_fluid_z.array());

		assert(!std::isnan(node_rhs_fluid_x[bucket_idx].sum()));
		assert(!std::isnan(node_rhs_fluid_y[bucket_idx].sum()));
		assert(!std::isnan(node_rhs_fluid_z[bucket_idx].sum()));
	});
}

void LinearizedImplicitEuler::addSolidDrag( TwoDScene& scene,
        const std::vector< VectorXs >& node_vel_x,
        const std::vector< VectorXs >& node_vel_y,
        const std::vector< VectorXs >& node_vel_z,
        std::vector< VectorXs >& node_fluid_vel_x,
        std::vector< VectorXs >& node_fluid_vel_y,
        std::vector< VectorXs >& node_fluid_vel_z )
{
	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			node_fluid_vel_x[bucket_idx][i] += m_node_hdvm_x[bucket_idx][i] * node_vel_x[bucket_idx][i] * m_node_inv_mfhdvm_x[bucket_idx][i];
			node_fluid_vel_y[bucket_idx][i] += m_node_hdvm_y[bucket_idx][i] * node_vel_y[bucket_idx][i] * m_node_inv_mfhdvm_y[bucket_idx][i];
			node_fluid_vel_z[bucket_idx][i] += m_node_hdvm_z[bucket_idx][i] * node_vel_z[bucket_idx][i] * m_node_inv_mfhdvm_z[bucket_idx][i];
		}
	});
}

void LinearizedImplicitEuler::addSolidDragRHS( TwoDScene& scene,
        const std::vector< VectorXs >& node_vel_x,
        const std::vector< VectorXs >& node_vel_y,
        const std::vector< VectorXs >& node_vel_z,
        std::vector< VectorXs >& node_rhs_x,
        std::vector< VectorXs >& node_rhs_y,
        std::vector< VectorXs >& node_rhs_z )
{
	const Sorter& buckets = scene.getParticleBuckets();

	const std::vector< VectorXs >& node_mass_x = scene.getNodeMassX();
	const std::vector< VectorXs >& node_mass_y = scene.getNodeMassY();
	const std::vector< VectorXs >& node_mass_z = scene.getNodeMassZ();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			node_rhs_x[bucket_idx][i] += m_node_mshdvm_hdvm_x[bucket_idx][i] * node_mass_x[bucket_idx][i] * node_vel_x[bucket_idx][i];
			node_rhs_y[bucket_idx][i] += m_node_mshdvm_hdvm_y[bucket_idx][i] * node_mass_y[bucket_idx][i] * node_vel_y[bucket_idx][i];
			node_rhs_z[bucket_idx][i] += m_node_mshdvm_hdvm_z[bucket_idx][i] * node_mass_z[bucket_idx][i] * node_vel_z[bucket_idx][i];
		}
	});
}


void LinearizedImplicitEuler::addFluidDragRHS( TwoDScene& scene,
        const std::vector< VectorXs >& node_fluid_vel_x,
        const std::vector< VectorXs >& node_fluid_vel_y,
        const std::vector< VectorXs >& node_fluid_vel_z,
        std::vector< VectorXs >& node_rhs_x,
        std::vector< VectorXs >& node_rhs_y,
        std::vector< VectorXs >& node_rhs_z )
{
	const Sorter& buckets = scene.getParticleBuckets();

	const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
	const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
	const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			node_rhs_x[bucket_idx][i] += m_node_mfhdvm_hdvm_x[bucket_idx][i] * node_mass_fluid_x[bucket_idx][i] * node_fluid_vel_x[bucket_idx][i];
			node_rhs_y[bucket_idx][i] += m_node_mfhdvm_hdvm_y[bucket_idx][i] * node_mass_fluid_y[bucket_idx][i] * node_fluid_vel_y[bucket_idx][i];
			node_rhs_z[bucket_idx][i] += m_node_mfhdvm_hdvm_z[bucket_idx][i] * node_mass_fluid_z[bucket_idx][i] * node_fluid_vel_z[bucket_idx][i];
		}
	});
}

void LinearizedImplicitEuler::constructHessianPreProcess( TwoDScene& scene, const scalar& dt )
{
	scene.accumulateddUdxdx(m_triA, dt, 0);

	m_triA.erase(std::remove_if(m_triA.begin(), m_triA.end(), [] (const auto & info) {return info.value() == 0.0;}), m_triA.end());
}

void LinearizedImplicitEuler::constructHessianPostProcess( TwoDScene& scene, const scalar& dt)
{
	const int num_soft_elasto = scene.getNumSoftElastoParticles();

	tbb::parallel_sort(m_triA.begin(), m_triA.end(), [] (const Triplets & x, const Triplets & y) {
		return x.row() < y.row();
	});

	if ((int) m_triA_sup.size() != num_soft_elasto * 4) m_triA_sup.resize(num_soft_elasto * 4);

	memset(&m_triA_sup[0], 0, num_soft_elasto * 4 * sizeof(std::pair<int, int>));

	const int num_tris = m_triA.size();

	threadutils::for_each(0, num_tris, [&] (int pidx) {
		int G_ID = pidx;
		int G_ID_PREV = G_ID - 1;
		int G_ID_NEXT = G_ID + 1;

		unsigned int cell = m_triA[G_ID].row();
		unsigned int cell_prev = G_ID_PREV < 0 ? -1U : m_triA[G_ID_PREV].row();
		unsigned int cell_next = G_ID_NEXT >= num_tris ? -1U : m_triA[G_ID_NEXT].row();

		if (cell != cell_prev) {
			m_triA_sup[cell].first = G_ID;
		}

		if (cell != cell_next) {
			m_triA_sup[cell].second = G_ID_NEXT;
		}
	});

	if (scene.getLiquidInfo().use_group_precondition) {
		prepareGroupPrecondition(scene, m_node_Cs_x, m_node_Cs_y, m_node_Cs_z, dt);
	}


}

void LinearizedImplicitEuler::constructAngularHessianPreProcess( TwoDScene& scene, const scalar& dt )
{
	scene.accumulateAngularddUdxdx(m_angular_triA, dt, 0);

	m_angular_triA.erase(std::remove_if(m_angular_triA.begin(), m_angular_triA.end(), [] (const auto & info) {return info.value() == 0.0;}), m_angular_triA.end());
}

void LinearizedImplicitEuler::constructAngularHessianPostProcess( TwoDScene& scene, const scalar& )
{
	const int num_soft_elasto = scene.getNumSoftElastoParticles();

	tbb::parallel_sort(m_angular_triA.begin(), m_angular_triA.end(), [] (const Triplets & x, const Triplets & y) {
		return x.row() < y.row();
	});

	if ((int) m_angular_triA_sup.size() != num_soft_elasto) m_angular_triA_sup.resize(num_soft_elasto);

	memset(&m_angular_triA_sup[0], 0, num_soft_elasto * sizeof(std::pair<int, int>));

	const int num_tris = m_angular_triA.size();

	threadutils::for_each(0, num_tris, [&] (int pidx) {
		int G_ID = pidx;
		int G_ID_PREV = G_ID - 1;
		int G_ID_NEXT = G_ID + 1;

		unsigned int cell = m_angular_triA[G_ID].row();
		unsigned int cell_prev = G_ID_PREV < 0 ? -1U : m_angular_triA[G_ID_PREV].row();
		unsigned int cell_next = G_ID_NEXT >= num_tris ? -1U : m_angular_triA[G_ID_NEXT].row();

		if (cell != cell_prev) {
			m_angular_triA_sup[cell].first = G_ID;
		}

		if (cell != cell_next) {
			m_angular_triA_sup[cell].second = G_ID_NEXT;
		}
	});
}

void LinearizedImplicitEuler::constructHDV( TwoDScene& scene, const scalar& dt )
{
	const std::vector< VectorXs >& node_mass_x = scene.getNodeMassX();
	const std::vector< VectorXs >& node_mass_y = scene.getNodeMassY();
	const std::vector< VectorXs >& node_mass_z = scene.getNodeMassZ();

	const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
	const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
	const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

	const std::vector< VectorXs >& node_vel_x = scene.getNodeVelocityX();
	const std::vector< VectorXs >& node_vel_y = scene.getNodeVelocityY();
	const std::vector< VectorXs >& node_vel_z = scene.getNodeVelocityZ();

	const std::vector< VectorXs >& node_vel_fluid_x = scene.getNodeFluidVelocityX();
	const std::vector< VectorXs >& node_vel_fluid_y = scene.getNodeFluidVelocityY();
	const std::vector< VectorXs >& node_vel_fluid_z = scene.getNodeFluidVelocityZ();

	const std::vector< VectorXs >& node_vol_fluid_x = scene.getNodeFluidVolX();
	const std::vector< VectorXs >& node_vol_fluid_y = scene.getNodeFluidVolY();
	const std::vector< VectorXs >& node_vol_fluid_z = scene.getNodeFluidVolZ();

	const std::vector< VectorXs >& node_vol_x = scene.getNodeVolX();
	const std::vector< VectorXs >& node_vol_y = scene.getNodeVolY();
	const std::vector< VectorXs >& node_vol_z = scene.getNodeVolZ();

	const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
	const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
	const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();

	const std::vector< VectorXs >& node_sat_x = scene.getNodeSaturationX();
	const std::vector< VectorXs >& node_sat_y = scene.getNodeSaturationY();
	const std::vector< VectorXs >& node_sat_z = scene.getNodeSaturationZ();

	const std::vector< VectorXs >& orient_x = scene.getNodeOrientationX();
	const std::vector< VectorXs >& orient_y = scene.getNodeOrientationY();
	const std::vector< VectorXs >& orient_z = scene.getNodeOrientationZ();

	const std::vector< VectorXs >& shape_factor_x = scene.getNodeShapeFactorX();
	const std::vector< VectorXs >& shape_factor_y = scene.getNodeShapeFactorY();
	const std::vector< VectorXs >& shape_factor_z = scene.getNodeShapeFactorZ();

	const Sorter& buckets = scene.getParticleBuckets();

	if (scene.getLiquidInfo().drag_by_air)
	{
		buckets.for_each_bucket([&] (int bucket_idx) {
			if (!scene.isBucketActivated(bucket_idx)) return;

			const int num_nodes = scene.getNumNodes(bucket_idx);

			for (int i = 0; i < num_nodes; ++i)
			{
				scalar dc = scene.getDragCoeffWithOrientation(node_psi_x[bucket_idx][i], 1.0, node_vel_x[bucket_idx][i], orient_x[bucket_idx].segment<3>(i * 3), shape_factor_x[bucket_idx][i], 0, 1);
				scalar V_m = node_vol_x[bucket_idx][i];
				scalar hdV = dt * dc * V_m;

				m_node_damped_x[bucket_idx][i] = node_mass_x[bucket_idx][i] + hdV;

				dc = scene.getDragCoeffWithOrientation(node_psi_y[bucket_idx][i], 1.0, node_vel_y[bucket_idx][i], orient_y[bucket_idx].segment<3>(i * 3), shape_factor_y[bucket_idx][i], 1, 1);
				V_m = node_vol_y[bucket_idx][i];
				hdV = dt * dc * V_m;

				m_node_damped_y[bucket_idx][i] = node_mass_y[bucket_idx][i] + hdV;

				dc = scene.getDragCoeffWithOrientation(node_psi_z[bucket_idx][i], 1.0, node_vel_z[bucket_idx][i], orient_z[bucket_idx].segment<3>(i * 3), shape_factor_z[bucket_idx][i], 2, 1);
				V_m = node_vol_z[bucket_idx][i];
				hdV = dt * dc * V_m;

				m_node_damped_z[bucket_idx][i] = node_mass_z[bucket_idx][i] + hdV;
			}
		});
	} else {
		m_node_damped_x = node_mass_x;
		m_node_damped_y = node_mass_y;
		m_node_damped_z = node_mass_z;
	}

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			scalar sat = node_sat_x[bucket_idx][i];
			scalar dc = scene.getDragCoeffWithOrientation(node_psi_x[bucket_idx][i], sat, node_vel_fluid_x[bucket_idx][i] - node_vel_x[bucket_idx][i], orient_x[bucket_idx].segment<3>(i * 3), shape_factor_x[bucket_idx][i], 0, 0);
			scalar V_m = node_vol_x[bucket_idx][i] + node_vol_fluid_x[bucket_idx][i];
			scalar hdV = dt * dc * V_m;

			m_node_hdvm_x[bucket_idx][i] = hdV;

			scalar mfhdvm = node_mass_fluid_x[bucket_idx][i] + hdV;
			scalar mshdvm = m_node_damped_x[bucket_idx][i] + hdV;

			m_node_inv_mfhdvm_x[bucket_idx][i] = (mfhdvm > 1e-20) ? (1.0 / mfhdvm) : 1.0;
			m_node_mfhdvm_hdvm_x[bucket_idx][i] = (mfhdvm > 1e-20) ? (hdV / mfhdvm) : 1.0;
			m_node_mshdvm_hdvm_x[bucket_idx][i] = (mshdvm > 1e-20) ? (hdV / mshdvm) : 1.0;

			sat = node_sat_y[bucket_idx][i];
			dc = scene.getDragCoeffWithOrientation(node_psi_y[bucket_idx][i], sat, node_vel_fluid_y[bucket_idx][i] - node_vel_y[bucket_idx][i], orient_y[bucket_idx].segment<3>(i * 3), shape_factor_y[bucket_idx][i], 1, 0);
			V_m = node_vol_y[bucket_idx][i] + node_vol_fluid_y[bucket_idx][i];
			hdV = dt * dc * V_m;

			m_node_hdvm_y[bucket_idx][i] = hdV;

			mfhdvm = node_mass_fluid_y[bucket_idx][i] + hdV;
			mshdvm = m_node_damped_y[bucket_idx][i] + hdV;

			m_node_inv_mfhdvm_y[bucket_idx][i] = (mfhdvm > 1e-20) ? (1.0 / mfhdvm) : 1.0;
			m_node_mfhdvm_hdvm_y[bucket_idx][i] = (mfhdvm > 1e-20) ? (hdV / mfhdvm) : 1.0;
			m_node_mshdvm_hdvm_y[bucket_idx][i] = (mshdvm > 1e-20) ? (hdV / mshdvm) : 1.0;

			sat = node_sat_z[bucket_idx][i];
			dc = scene.getDragCoeffWithOrientation(node_psi_z[bucket_idx][i], sat, node_vel_fluid_z[bucket_idx][i] - node_vel_z[bucket_idx][i], orient_z[bucket_idx].segment<3>(i * 3), shape_factor_z[bucket_idx][i], 2, 0);
			V_m = node_vol_z[bucket_idx][i] + node_vol_fluid_z[bucket_idx][i];
			hdV = dt * dc * V_m;

			m_node_hdvm_z[bucket_idx][i] = hdV;

			mfhdvm = node_mass_fluid_z[bucket_idx][i] + hdV;
			mshdvm = m_node_damped_z[bucket_idx][i] + hdV;

			m_node_inv_mfhdvm_z[bucket_idx][i] = (mfhdvm > 1e-20) ? (1.0 / mfhdvm) : 1.0;
			m_node_mfhdvm_hdvm_z[bucket_idx][i] = (mfhdvm > 1e-20) ? (hdV / mfhdvm) : 1.0;
			m_node_mshdvm_hdvm_z[bucket_idx][i] = (mshdvm > 1e-20) ? (hdV / mshdvm) : 1.0;
		}
	});
}

void LinearizedImplicitEuler::constructMsDVs( TwoDScene& scene )
{
	const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
	const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
	const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			scalar P = m_node_mfhdvm_hdvm_x[bucket_idx][i];
			scalar Cs = m_node_damped_x[bucket_idx][i] + P * node_mass_fluid_x[bucket_idx][i];
			m_node_Cs_x[bucket_idx][i] = Cs;
			if (Cs > 1e-20) {
				m_node_inv_Cs_x[bucket_idx][i] = 1.0 / Cs;
			} else {
				m_node_inv_Cs_x[bucket_idx][i] = 1.0;
			}

			P = m_node_mfhdvm_hdvm_y[bucket_idx][i];
			Cs = m_node_damped_y[bucket_idx][i] + P * node_mass_fluid_y[bucket_idx][i];
			m_node_Cs_y[bucket_idx][i] = Cs;
			if (Cs > 1e-20) {
				m_node_inv_Cs_y[bucket_idx][i] = 1.0 / Cs;
			} else {
				m_node_inv_Cs_y[bucket_idx][i] = 1.0;
			}

			P = m_node_mfhdvm_hdvm_z[bucket_idx][i];
			Cs = m_node_damped_z[bucket_idx][i] + P * node_mass_fluid_z[bucket_idx][i];
			m_node_Cs_z[bucket_idx][i] = Cs;
			if (Cs > 1e-20) {
				m_node_inv_Cs_z[bucket_idx][i] = 1.0 / Cs;
			} else {
				m_node_inv_Cs_z[bucket_idx][i] = 1.0;
			}
		}
	});

//	mathutils::print_histogram_analysis(m_node_Cs_x, 10, "Cs_x", true);
//	mathutils::print_histogram_analysis(m_node_Cs_y, 10, "Cs_y", true);
//	mathutils::print_histogram_analysis(m_node_Cs_z, 10, "Cs_z", true);
}

void LinearizedImplicitEuler::constructPsiSF( TwoDScene& scene )
{
	const std::vector< VectorXs >& node_psi_x = scene.getNodePsiX();
	const std::vector< VectorXs >& node_psi_y = scene.getNodePsiY();
	const std::vector< VectorXs >& node_psi_z = scene.getNodePsiZ();

	const std::vector< VectorXs >& node_mass_x = scene.getNodeMassX();
	const std::vector< VectorXs >& node_mass_y = scene.getNodeMassY();
	const std::vector< VectorXs >& node_mass_z = scene.getNodeMassZ();

	const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
	const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
	const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			m_node_psi_fs_x[bucket_idx][i] = ((1.0 - node_psi_x[bucket_idx][i]) * m_node_inv_C_x[bucket_idx][i] +
			                                  node_psi_x[bucket_idx][i] * m_node_inv_Cs_x[bucket_idx][i] * m_node_mfhdvm_hdvm_x[bucket_idx][i]) * node_mass_fluid_x[bucket_idx][i];
			m_node_psi_sf_x[bucket_idx][i] = (node_psi_x[bucket_idx][i] * m_node_inv_Cs_x[bucket_idx][i] +
			                                  (1.0 - node_psi_x[bucket_idx][i]) * m_node_inv_C_x[bucket_idx][i] * m_node_mshdvm_hdvm_x[bucket_idx][i]) * node_mass_x[bucket_idx][i];

			m_node_psi_fs_y[bucket_idx][i] = ((1.0 - node_psi_y[bucket_idx][i]) * m_node_inv_C_y[bucket_idx][i] +
			                                  node_psi_y[bucket_idx][i] * m_node_inv_Cs_y[bucket_idx][i] * m_node_mfhdvm_hdvm_y[bucket_idx][i]) * node_mass_fluid_y[bucket_idx][i];
			m_node_psi_sf_y[bucket_idx][i] = (node_psi_y[bucket_idx][i] * m_node_inv_Cs_y[bucket_idx][i] +
			                                  (1.0 - node_psi_y[bucket_idx][i]) * m_node_inv_C_y[bucket_idx][i] * m_node_mshdvm_hdvm_y[bucket_idx][i]) * node_mass_y[bucket_idx][i];

			m_node_psi_fs_z[bucket_idx][i] = ((1.0 - node_psi_z[bucket_idx][i]) * m_node_inv_C_z[bucket_idx][i] +
			                                  node_psi_z[bucket_idx][i] * m_node_inv_Cs_z[bucket_idx][i] * m_node_mfhdvm_hdvm_z[bucket_idx][i]) * node_mass_fluid_z[bucket_idx][i];
			m_node_psi_sf_z[bucket_idx][i] = (node_psi_z[bucket_idx][i] * m_node_inv_Cs_z[bucket_idx][i] +
			                                  (1.0 - node_psi_z[bucket_idx][i]) * m_node_inv_C_z[bucket_idx][i] * m_node_mshdvm_hdvm_z[bucket_idx][i]) * node_mass_z[bucket_idx][i];
		}
	});

}

void LinearizedImplicitEuler::constructInvMDV( TwoDScene& scene )
{
	const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
	const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
	const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

	const Sorter& buckets = scene.getParticleBuckets();

	buckets.for_each_bucket([&] (int bucket_idx) {
		if (!scene.isBucketActivated(bucket_idx)) return;

		const int num_nodes = scene.getNumNodes(bucket_idx);

		for (int i = 0; i < num_nodes; ++i)
		{
			scalar Q = m_node_mshdvm_hdvm_x[bucket_idx][i];
			scalar C = node_mass_fluid_x[bucket_idx][i] + Q * m_node_damped_x[bucket_idx][i];
			if (C > 1e-20) {
				m_node_inv_C_x[bucket_idx][i] = 1.0 / C;
			} else {
				m_node_inv_C_x[bucket_idx][i] = 1.0;
			}

			Q = m_node_mshdvm_hdvm_y[bucket_idx][i];
			C = node_mass_fluid_y[bucket_idx][i] + Q * m_node_damped_y[bucket_idx][i];
			if (C > 1e-20) {
				m_node_inv_C_y[bucket_idx][i] = 1.0 / C;
			} else {
				m_node_inv_C_y[bucket_idx][i] = 1.0;
			}

			Q = m_node_mshdvm_hdvm_z[bucket_idx][i];
			C = node_mass_fluid_z[bucket_idx][i] + Q * m_node_damped_z[bucket_idx][i];
			if (C > 1e-20) {
				m_node_inv_C_z[bucket_idx][i] = 1.0 / C;
			} else {
				m_node_inv_C_z[bucket_idx][i] = 1.0;
			}
		}
	});
}

scalar LinearizedImplicitEuler::computeDivergence( TwoDScene& scene )
{
	std::vector< VectorXs > node_ic;
	allocateCenterNodeVectors(scene, node_ic);

	std::vector< scalar > bucket_rms(node_ic.size());

	pressure::constructNodeIncompressibleCondition(scene, node_ic, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z, m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z);

	scene.getParticleBuckets().for_each_bucket([&] (int bucket_idx) {
		bucket_rms[bucket_idx] = node_ic[bucket_idx].squaredNorm();
	});

	scalar rms = 0.0;
	int count = 0;
	const int num_buckets = node_ic.size();
	for (int i = 0; i < num_buckets; ++i) {
		rms += bucket_rms[i];
		count += (int) node_ic[i].size();
	}

	count = std::max(1, count);

	return sqrt( rms / (scalar) count );
}

bool LinearizedImplicitEuler::stepImplicitElastoDiagonalPCR( TwoDScene& scene, scalar dt )
{
	int ndof_elasto = scene.getNumSoftElastoParticles() * 4;
	const Sorter& buckets = scene.getParticleBuckets();

	if (ndof_elasto == 0) return true;

	scalar res_norm_0 = lengthNodeVectors(m_node_rhs_x, m_node_rhs_y, m_node_rhs_z);
	scalar res_norm_1 = m_angular_moment_buffer.norm();

	if (res_norm_0 > m_pcg_criterion) {
		constructHessianPreProcess(scene, dt);
		constructHessianPostProcess(scene, dt);

		allocateNodeVectors(scene, m_node_r_x, m_node_r_y, m_node_r_z);
		allocateNodeVectors(scene, m_node_z_x, m_node_z_y, m_node_z_z);
		allocateNodeVectors(scene, m_node_p_x, m_node_p_y, m_node_p_z);
		allocateNodeVectors(scene, m_node_q_x, m_node_q_y, m_node_q_z);
		allocateNodeVectors(scene, m_node_w_x, m_node_w_y, m_node_w_z);
		allocateNodeVectors(scene, m_node_t_x, m_node_t_y, m_node_t_z);

		// initial residual = b - Ax
		performGlobalMultiply(scene, dt,
		                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
		                      m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z,
		                      m_node_z_x, m_node_z_y, m_node_z_z);

		buckets.for_each_bucket([&] (int bucket_idx) {
			m_node_z_x[bucket_idx] = m_node_rhs_x[bucket_idx] - m_node_z_x[bucket_idx];
			m_node_z_y[bucket_idx] = m_node_rhs_y[bucket_idx] - m_node_z_y[bucket_idx];
			m_node_z_z[bucket_idx] = m_node_rhs_z[bucket_idx] - m_node_z_z[bucket_idx];
		});

		scalar res_norm = lengthNodeVectors(m_node_z_x, m_node_z_y, m_node_z_z) / res_norm_0;

		int iter = 0;

		if (res_norm < m_pcg_criterion) {
			std::cout << "[pcr total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_0) << "/" << (m_pcg_criterion * res_norm_0)
			          << "]" << std::endl;
		} else {
			// Solve Mr=z
			if (scene.getLiquidInfo().use_group_precondition) {
				performGroupedLocalSolve(scene, m_node_z_x, m_node_z_y, m_node_z_z,
				                         m_node_r_x, m_node_r_y, m_node_r_z);
			} else {
				performInvLocalSolve(scene, m_node_z_x, m_node_z_y, m_node_z_z,
				                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
				                     m_node_r_x, m_node_r_y, m_node_r_z);
			}
			// p = r
			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_p_x[bucket_idx] = m_node_r_x[bucket_idx];
				m_node_p_y[bucket_idx] = m_node_r_y[bucket_idx];
				m_node_p_z[bucket_idx] = m_node_r_z[bucket_idx];
			});

			// t = z
			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_t_x[bucket_idx] = m_node_z_x[bucket_idx];
				m_node_t_y[bucket_idx] = m_node_z_y[bucket_idx];
				m_node_t_z[bucket_idx] = m_node_z_z[bucket_idx];
			});

			// w = Ar
			performGlobalMultiply(scene, dt,
			                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
			                      m_node_r_x, m_node_r_y, m_node_r_z,
			                      m_node_w_x, m_node_w_y, m_node_w_z);

			// rho = (r, w)
			scalar rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_w_x, m_node_w_y, m_node_w_z);

			// q = Ap
			performGlobalMultiply(scene, dt,
			                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
			                      m_node_p_x, m_node_p_y, m_node_p_z,
			                      m_node_q_x, m_node_q_y, m_node_q_z);

			if (scene.getLiquidInfo().use_group_precondition) {
				performGroupedLocalSolve(scene, m_node_q_x, m_node_q_y, m_node_q_z,
				                         m_node_z_x, m_node_z_y, m_node_z_z);
			} else {
				// Mz=q
				performInvLocalSolve(scene, m_node_q_x, m_node_q_y, m_node_q_z,
				                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
				                     m_node_z_x, m_node_z_y, m_node_z_z);
			}


			// alpha = rho / (q, z)
			scalar alpha = rho / dotNodeVectors(m_node_q_x, m_node_q_y, m_node_q_z, m_node_z_x, m_node_z_y, m_node_z_z);

			// x = x + alpha * p
			// r = r - alpha * z
			// t = t - alpha * q
			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_v_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
				m_node_r_x[bucket_idx] -= m_node_z_x[bucket_idx] * alpha;
				m_node_t_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
				m_node_v_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
				m_node_r_y[bucket_idx] -= m_node_z_y[bucket_idx] * alpha;
				m_node_t_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
				m_node_v_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
				m_node_r_z[bucket_idx] -= m_node_z_z[bucket_idx] * alpha;
				m_node_t_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
			});

			res_norm = lengthNodeVectors(m_node_t_x, m_node_t_y, m_node_t_z) / res_norm_0;

			const scalar rho_criterion = (m_pcg_criterion * res_norm_0) * (m_pcg_criterion * res_norm_0);

			scalar rho_old, beta;
			for (; iter < m_maxiters && res_norm > m_pcg_criterion && rho > rho_criterion; ++iter)
			{
				rho_old = rho;

				// w = Ar
				performGlobalMultiply(scene, dt,
				                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
				                      m_node_r_x, m_node_r_y, m_node_r_z,
				                      m_node_w_x, m_node_w_y, m_node_w_z);

				// rho = (r, w)
				rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_w_x, m_node_w_y, m_node_w_z);

				beta = rho / rho_old;

				// p = beta * p + r
				// q = beta * q + w
				buckets.for_each_bucket([&] (int bucket_idx) {
					m_node_p_x[bucket_idx] = m_node_r_x[bucket_idx] + m_node_p_x[bucket_idx] * beta;
					m_node_p_y[bucket_idx] = m_node_r_y[bucket_idx] + m_node_p_y[bucket_idx] * beta;
					m_node_p_z[bucket_idx] = m_node_r_z[bucket_idx] + m_node_p_z[bucket_idx] * beta;
					m_node_q_x[bucket_idx] = m_node_w_x[bucket_idx] + m_node_q_x[bucket_idx] * beta;
					m_node_q_y[bucket_idx] = m_node_w_y[bucket_idx] + m_node_q_y[bucket_idx] * beta;
					m_node_q_z[bucket_idx] = m_node_w_z[bucket_idx] + m_node_q_z[bucket_idx] * beta;
				});

				// Mz = q
				if (scene.getLiquidInfo().use_group_precondition) {
					performGroupedLocalSolve(scene, m_node_q_x, m_node_q_y, m_node_q_z,
					                         m_node_z_x, m_node_z_y, m_node_z_z);
				} else {
					// Mz=q
					performInvLocalSolve(scene, m_node_q_x, m_node_q_y, m_node_q_z,
					                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
					                     m_node_z_x, m_node_z_y, m_node_z_z);
				}

				// alpha = rho / (q, z)
				alpha = rho / dotNodeVectors(m_node_q_x, m_node_q_y, m_node_q_z, m_node_z_x, m_node_z_y, m_node_z_z);

				// x = x + alpha * p
				// r = r - alpha * z
				// t = t - alpha * q
				buckets.for_each_bucket([&] (int bucket_idx) {
					m_node_v_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
					m_node_r_x[bucket_idx] -= m_node_z_x[bucket_idx] * alpha;
					m_node_t_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
					m_node_v_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
					m_node_r_y[bucket_idx] -= m_node_z_y[bucket_idx] * alpha;
					m_node_t_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
					m_node_v_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
					m_node_r_z[bucket_idx] -= m_node_z_z[bucket_idx] * alpha;
					m_node_t_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
				});

				res_norm = lengthNodeVectors(m_node_t_x, m_node_t_y, m_node_t_z) / res_norm_0;
				if (scene.getLiquidInfo().iteration_print_step > 0 && iter % scene.getLiquidInfo().iteration_print_step == 0)
					std::cout << "[pcr total iter: " << iter
					          << ", res: " << res_norm << "/" << m_pcg_criterion
					          << ", abs. res: " << (res_norm * res_norm_0) << "/" << (m_pcg_criterion * res_norm_0)
					          << ", rho: " << (rho / (res_norm_0 * res_norm_0)) << "/" << (rho_criterion / (res_norm_0 * res_norm_0))
					          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;
			}

			std::cout << "[pcr total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_0) << "/" << (m_pcg_criterion * res_norm_0)
			          << ", rho: " << (rho / (res_norm_0 * res_norm_0)) << "/" << (rho_criterion / (res_norm_0 * res_norm_0))
			          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;
		}
	}

	if (res_norm_1 > m_pcg_criterion) {
		const int num_elasto = scene.getNumSoftElastoParticles();

		constructAngularHessianPreProcess(scene, dt);
		constructAngularHessianPostProcess(scene, dt);

		m_angular_r.resize(num_elasto);
		m_angular_z.resize(num_elasto);
		m_angular_p.resize(num_elasto);
		m_angular_q.resize(num_elasto);
		m_angular_w.resize(num_elasto);
		m_angular_t.resize(num_elasto);

		m_angular_r.setZero();
		m_angular_z.setZero();
		m_angular_p.setZero();
		m_angular_q.setZero();
		m_angular_w.setZero();
		m_angular_t.setZero();

		performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_v_plus_buffer, m_angular_z);

		m_angular_z = m_angular_moment_buffer - m_angular_z;

		scalar res_norm = m_angular_z.norm() / res_norm_1;

		int iter = 0;

		if (res_norm < m_pcg_criterion ) {
			std::cout << "[angular pcr total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_1) << "/" << (m_pcg_criterion * res_norm_1)
			          << "]" << std::endl;
		} else {

			// Solve Mr=z
			performLocalSolveTwist(scene, m_angular_z, scene.getM(), m_angular_r);

			// p = r
			m_angular_p = m_angular_r;

			// t = z
			m_angular_t = m_angular_z;

			// w = Ar
			performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_r, m_angular_w);

			// rho = (r, w)
			scalar rho = m_angular_r.dot(m_angular_w);

			// q = Ap
			performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_p, m_angular_q);

			// Mz=q
			performLocalSolveTwist(scene, m_angular_q, scene.getM(), m_angular_z);

			// alpha = rho / (q, z)
			scalar alpha = rho / m_angular_q.dot(m_angular_z);

			// x = x + alpha * p
			// r = r - alpha * z
			// t = t - alpha * q
			m_angular_v_plus_buffer += m_angular_p * alpha;
			m_angular_r -= m_angular_z * alpha;
			m_angular_t -= m_angular_q * alpha;

			res_norm = m_angular_t.norm() / res_norm_1;

			const scalar rho_criterion = (m_pcg_criterion * res_norm_1) * (m_pcg_criterion * res_norm_1);

			scalar rho_old, beta;
			for (; iter < m_maxiters && res_norm > m_pcg_criterion && rho > rho_criterion; ++iter)
			{
				rho_old = rho;

				// w = Ar
				performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_r, m_angular_w);

				// rho = (r, w)
				rho = m_angular_r.dot(m_angular_w);

				beta = rho / rho_old;

				// p = beta * p + r
				// q = beta * q + w
				m_angular_p = m_angular_r + m_angular_p * beta;
				m_angular_q = m_angular_w + m_angular_q * beta;

				// Mz = q
				performLocalSolveTwist(scene, m_angular_q, scene.getM(), m_angular_z);

				// alpha = rho / (q, z)
				alpha = rho / m_angular_q.dot(m_angular_z);

				// x = x + alpha * p
				// r = r - alpha * z
				// t = t - alpha * q
				m_angular_v_plus_buffer += m_angular_p * alpha;
				m_angular_r -= m_angular_z * alpha;
				m_angular_t -= m_angular_q * alpha;

				res_norm = m_angular_t.norm() / res_norm_1;

				if (scene.getLiquidInfo().iteration_print_step > 0 && iter % scene.getLiquidInfo().iteration_print_step == 0)
					std::cout << "[angular pcr total iter: " << iter
					          << ", res: " << res_norm << "/" << m_pcg_criterion
					          << ", abs. res: " << (res_norm * res_norm_1) << "/" << (m_pcg_criterion * res_norm_1)
					          << ", rho: " << (rho / (res_norm_1 * res_norm_1)) << "/" << (rho_criterion / (res_norm_1 * res_norm_1))
					          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;

			}

			std::cout << "[angular pcr total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_1) << "/" << (m_pcg_criterion * res_norm_1)
			          << ", rho: " << (rho / (res_norm_1 * res_norm_1)) << "/" << (rho_criterion / (res_norm_1 * res_norm_1))
			          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;
		}
	}

	return true;
}

void LinearizedImplicitEuler::performGroupedLocalSolve( const TwoDScene& scene,
        const std::vector< VectorXs >& node_rhs_x,
        const std::vector< VectorXs >& node_rhs_y,
        const std::vector< VectorXs >& node_rhs_z,
        std::vector< VectorXs >& out_node_vec_x,
        std::vector< VectorXs >& out_node_vec_y,
        std::vector< VectorXs >& out_node_vec_z )
{
	const std::vector< VectorXi >& groups = scene.getSolveGroup();

	const int num_groups = (int) groups.size();

	VectorXs rhs_buffer(scene.getNumSoftElastoParticles() * 4);
	rhs_buffer.setZero();

	VectorXs sol_buffer(scene.getNumSoftElastoParticles() * 4);
	sol_buffer.setZero();

	mapNodeToSoftParticles( scene, node_rhs_x, node_rhs_y, node_rhs_z, rhs_buffer );

	threadutils::for_each(0, num_groups, [&] (int igroup) {
		const VectorXi& members = groups[igroup];
		const int num_members = members.size();

		VectorXs group_rhs(num_members * 3);
		for (int i = 0; i < num_members; ++i)
		{
			group_rhs.segment<3>(i * 3) = rhs_buffer.segment<3>(members[i] * 4);
		}

		VectorXs group_sol = m_group_preconditioners[igroup]->solve(group_rhs);

		for (int i = 0; i < num_members; ++i)
		{
			sol_buffer.segment<3>(members[i] * 4) = group_sol.segment<3>(i * 3);
		}
	});

	mapSoftParticlesToNode( scene, out_node_vec_x, out_node_vec_y, out_node_vec_z, sol_buffer);
}

void LinearizedImplicitEuler::prepareGroupPrecondition(
    const TwoDScene& scene,
    const std::vector< VectorXs >& node_m_x,
    const std::vector< VectorXs >& node_m_y,
    const std::vector< VectorXs >& node_m_z,
    const scalar& dt )
{
	const std::vector< VectorXi >& groups = scene.getSolveGroup();

	const int num_groups = (int) groups.size();

	m_group_preconditioners.resize(num_groups);

	VectorXs mass_buffer(scene.getNumSoftElastoParticles() * 4);
	mass_buffer.setZero();
	// map drag + mass vector back to particles
	mapNodeToSoftParticles( scene, node_m_x, node_m_y, node_m_z, mass_buffer );

	threadutils::for_each(0, num_groups, [&] (int igroup) {
		const VectorXi& members = groups[igroup];
		std::unordered_map<int, int> finder;

		const int num_members = members.size();

		for (int i = 0; i < num_members; ++i) {
			finder[members[i]] = i;
		}

		TripletXs tri_sub_A;

		for (auto p : finder) {
			for (int r = 0; r < 3; ++r) {
				int i = p.first * 4 + r;
				const int idata_start = m_triA_sup[i].first;
				const int idata_end = m_triA_sup[i].second;

				for (int j = idata_start; j < idata_end; ++j)
				{
					const Triplets& tri = m_triA[j];
					const int qidx = tri.col() / 4;
					const int s = tri.col() - qidx * 4;
					if (s >= 3) continue;

					auto q = finder.find(qidx);
					if (q == finder.end()) continue;

					tri_sub_A.push_back(Triplets(p.second * 3 + r, q->second * 3 + s, tri.value() * dt * dt));
				}
			}
		}

		for (int i = 0; i < num_members; ++i)
		{
			const int pidx = members[i];

			for (int r = 0; r < 3; ++r) {
				tri_sub_A.push_back(Triplets(i * 3 + r, i * 3 + r, mass_buffer[pidx * 4 + r]));
			}
		}

		SparseXs sub_A(num_members * 3, num_members * 3);
		sub_A.setFromTriplets(tri_sub_A.begin(), tri_sub_A.end());

		m_group_preconditioners[igroup] = std::make_shared< Eigen::SimplicialLDLT< SparseXs > >(sub_A);
	});
}

bool LinearizedImplicitEuler::stepImplicitElastoLagrangian( TwoDScene& scene, scalar dt )
{
	int ndof_elasto = scene.getNumSoftElastoParticles() * 4;
	if (ndof_elasto == 0) return true;

	scalar res_norm_0 = m_lagrangian_rhs.norm();

	if (res_norm_0 > m_pcg_criterion) {
		// build Hessian
		constructHessianPreProcess(scene, dt);
		constructHessianPostProcess(scene, dt);

		allocateLagrangianVectors(scene, m_r);
		allocateLagrangianVectors(scene, m_z);
		allocateLagrangianVectors(scene, m_p);
		allocateLagrangianVectors(scene, m_q);

		const VectorXs& m = scene.getM();

		performGlobalMultiply(scene, dt, m, m_v_plus, m_r);

		m_r = m_lagrangian_rhs - m_r;

		scalar res_norm = m_r.norm() / res_norm_0;

		int iter = 0;

		if (res_norm < m_pcg_criterion) {
			std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "]" << std::endl;
		} else {
			performLocalSolve(scene, m_r, m, m_z);

			m_p = m_z;

			performGlobalMultiply(scene, dt, m, m_p, m_q);

			scalar rho = m_r.dot(m_z);

			scalar alpha = rho / m_p.dot(m_q);

			m_v_plus += m_p * alpha;
			m_r -= m_q * alpha;

			res_norm = m_r.norm() / res_norm_0;

			scalar rho_old, beta;
			for (; iter < m_maxiters && res_norm > m_pcg_criterion && rho > m_pcg_criterion * res_norm_0; ++iter)
			{
				rho_old = rho;

				performLocalSolve(scene, m_r, m, m_z);

				rho = m_r.dot(m_z);

				beta = rho / rho_old;

				m_p = m_z + m_p * beta;

				performGlobalMultiply(scene, dt, m, m_p, m_q);

				alpha = rho / m_p.dot(m_q);

				m_v_plus += m_p * alpha;
				m_r -= m_q * alpha;

				res_norm = m_r.norm() / res_norm_0;

				if (scene.getLiquidInfo().iteration_print_step > 0 && iter % scene.getLiquidInfo().iteration_print_step == 0)
					std::cout << "[pcg iter: " << iter << ", res: " << res_norm << "]" << std::endl;

			}

			std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "]" << std::endl;
		}
	}

	scene.getV().segment(0, ndof_elasto) = m_v_plus;
	return true;
}

bool LinearizedImplicitEuler::stepImplicitElastoDiagonalPCRCoSolve( TwoDScene& scene, scalar dt )
{
	return true;
}

bool LinearizedImplicitEuler::stepImplicitElastoDiagonalPCGCoSolve( TwoDScene& scene, scalar dt )
{
	const int num_elasto = scene.getNumSoftElastoParticles();
	const int ndof_elasto = num_elasto * 4;
	const Sorter& buckets = scene.getParticleBuckets();

	if (ndof_elasto == 0) return true;

	scalar res_norm_0 = lengthNodeVectors(m_node_rhs_x, m_node_rhs_y, m_node_rhs_z, m_angular_moment_buffer);
#ifdef PCG_VERBOSE
	std::cout << "[pcg total res0: " << res_norm_0 << "]" << std::endl;
#endif
	if (res_norm_0 > m_pcg_criterion) {
		// build Hessian
		constructHessianPreProcess(scene, dt);
		constructHessianPostProcess(scene, dt);

		allocateNodeVectors(scene, m_node_r_x, m_node_r_y, m_node_r_z);
		allocateNodeVectors(scene, m_node_z_x, m_node_z_y, m_node_z_z);
		allocateNodeVectors(scene, m_node_p_x, m_node_p_y, m_node_p_z);
		allocateNodeVectors(scene, m_node_q_x, m_node_q_y, m_node_q_z);

		m_angular_r.resize(num_elasto);
		m_angular_z.resize(num_elasto);
		m_angular_p.resize(num_elasto);
		m_angular_q.resize(num_elasto);

		m_angular_r.setZero();
		m_angular_z.setZero();
		m_angular_p.setZero();
		m_angular_q.setZero();

		performGlobalMultiply(scene, dt,
		                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
		                      m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z,
		                      m_node_r_x, m_node_r_y, m_node_r_z,
		                      scene.getM(), m_angular_v_plus_buffer, m_angular_r);

		buckets.for_each_bucket([&] (int bucket_idx) {
			m_node_r_x[bucket_idx] = m_node_rhs_x[bucket_idx] - m_node_r_x[bucket_idx];
			m_node_r_y[bucket_idx] = m_node_rhs_y[bucket_idx] - m_node_r_y[bucket_idx];
			m_node_r_z[bucket_idx] = m_node_rhs_z[bucket_idx] - m_node_r_z[bucket_idx];
		});

		m_angular_r = m_angular_moment_buffer - m_angular_r;

		scalar res_norm = lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_angular_r) / res_norm_0;

#ifdef PCG_VERBOSE
		std::cout << "[pcg total res: " << res_norm << "]" << std::endl;
#endif

		int iter = 0;

		if (res_norm < m_pcg_criterion) {
			std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "]" << std::endl;
		} else {
			performInvLocalSolve(scene, m_node_r_x, m_node_r_y, m_node_r_z,
			                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
			                     m_node_z_x, m_node_z_y, m_node_z_z);

			performLocalSolveTwist(scene, m_angular_r, scene.getM(), m_angular_z);

			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_p_x[bucket_idx] = m_node_z_x[bucket_idx];
				m_node_p_y[bucket_idx] = m_node_z_y[bucket_idx];
				m_node_p_z[bucket_idx] = m_node_z_z[bucket_idx];
			});

			m_angular_p = m_angular_z;

			performGlobalMultiply(scene, dt,
			                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
			                      m_node_p_x, m_node_p_y, m_node_p_z,
			                      m_node_q_x, m_node_q_y, m_node_q_z,
			                      scene.getM(), m_angular_p, m_angular_q);

			scalar rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_z_x, m_node_z_y, m_node_z_z, m_angular_r, m_angular_z);

			scalar alpha = rho / dotNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z, m_node_q_x, m_node_q_y, m_node_q_z, m_angular_p, m_angular_q);

			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_v_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
				m_node_r_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
				m_node_v_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
				m_node_r_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
				m_node_v_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
				m_node_r_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
			});

			m_angular_v_plus_buffer += m_angular_p * alpha;
			m_angular_r -= m_angular_q * alpha;

			res_norm = lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_angular_r) / res_norm_0;

			scalar rho_old, beta;
			for (; iter < m_maxiters && res_norm > m_pcg_criterion && rho > m_pcg_criterion * res_norm_0; ++iter)
			{
				rho_old = rho;

				performInvLocalSolve(scene, m_node_r_x, m_node_r_y, m_node_r_z,
				                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
				                     m_node_z_x, m_node_z_y, m_node_z_z);

				performLocalSolveTwist(scene, m_angular_r, scene.getM(), m_angular_z);

				rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_z_x, m_node_z_y, m_node_z_z, m_angular_r, m_angular_z);

				beta = rho / rho_old;

				buckets.for_each_bucket([&] (int bucket_idx) {
					m_node_p_x[bucket_idx] = m_node_z_x[bucket_idx] + m_node_p_x[bucket_idx] * beta;
					m_node_p_y[bucket_idx] = m_node_z_y[bucket_idx] + m_node_p_y[bucket_idx] * beta;
					m_node_p_z[bucket_idx] = m_node_z_z[bucket_idx] + m_node_p_z[bucket_idx] * beta;
				});

				m_angular_p = m_angular_z + m_angular_p * beta;

				performGlobalMultiply(scene, dt,
				                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
				                      m_node_p_x, m_node_p_y, m_node_p_z,
				                      m_node_q_x, m_node_q_y, m_node_q_z,
				                      scene.getM(), m_angular_p, m_angular_q);

				alpha = rho / dotNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z, m_node_q_x, m_node_q_y, m_node_q_z, m_angular_p, m_angular_q);

				buckets.for_each_bucket([&] (int bucket_idx) {
					m_node_v_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
					m_node_r_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
					m_node_v_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
					m_node_r_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
					m_node_v_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
					m_node_r_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
				});

				m_angular_v_plus_buffer += m_angular_p * alpha;
				m_angular_r -= m_angular_q * alpha;

				res_norm = lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_angular_r) / res_norm_0;

				if (scene.getLiquidInfo().iteration_print_step > 0 && iter % scene.getLiquidInfo().iteration_print_step == 0)
					std::cout << "[pcg iter: " << iter << ", res: " << res_norm << "]" << std::endl;

			}

			std::cout << "[pcg total iter: " << iter << ", res: " << res_norm << "]" << std::endl;
		}
	}

	return true;
}

bool LinearizedImplicitEuler::stepImplicitElastoDiagonalPCG( TwoDScene& scene, scalar dt )
{
	int ndof_elasto = scene.getNumSoftElastoParticles() * 4;
	const Sorter& buckets = scene.getParticleBuckets();

	if (ndof_elasto == 0) return true;

	scalar res_norm_0 = lengthNodeVectors(m_node_rhs_x, m_node_rhs_y, m_node_rhs_z);
	scalar res_norm_1 = m_angular_moment_buffer.norm();

	if (res_norm_0 > m_pcg_criterion) {
		// build Hessian
		constructHessianPreProcess(scene, dt);
		constructHessianPostProcess(scene, dt);

		allocateNodeVectors(scene, m_node_r_x, m_node_r_y, m_node_r_z);
		allocateNodeVectors(scene, m_node_z_x, m_node_z_y, m_node_z_z);
		allocateNodeVectors(scene, m_node_p_x, m_node_p_y, m_node_p_z);
		allocateNodeVectors(scene, m_node_q_x, m_node_q_y, m_node_q_z);

		performGlobalMultiply(scene, dt,
		                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
		                      m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z,
		                      m_node_r_x, m_node_r_y, m_node_r_z);

		buckets.for_each_bucket([&] (int bucket_idx) {
			m_node_r_x[bucket_idx] = m_node_rhs_x[bucket_idx] - m_node_r_x[bucket_idx];
			m_node_r_y[bucket_idx] = m_node_rhs_y[bucket_idx] - m_node_r_y[bucket_idx];
			m_node_r_z[bucket_idx] = m_node_rhs_z[bucket_idx] - m_node_r_z[bucket_idx];
		});

		scalar res_norm = lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z) / res_norm_0;

		int iter = 0;

		if (res_norm < m_pcg_criterion) {
			std::cout << "[pcg total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_0) << "/" << (m_pcg_criterion * res_norm_0)
			          << "]" << std::endl;
		} else {
			performInvLocalSolve(scene, m_node_r_x, m_node_r_y, m_node_r_z,
			                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
			                     m_node_z_x, m_node_z_y, m_node_z_z);

			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_p_x[bucket_idx] = m_node_z_x[bucket_idx];
				m_node_p_y[bucket_idx] = m_node_z_y[bucket_idx];
				m_node_p_z[bucket_idx] = m_node_z_z[bucket_idx];
			});

			performGlobalMultiply(scene, dt,
			                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
			                      m_node_p_x, m_node_p_y, m_node_p_z,
			                      m_node_q_x, m_node_q_y, m_node_q_z);

			scalar rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_z_x, m_node_z_y, m_node_z_z);

			scalar alpha = rho / dotNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z, m_node_q_x, m_node_q_y, m_node_q_z);

			buckets.for_each_bucket([&] (int bucket_idx) {
				m_node_v_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
				m_node_r_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
				m_node_v_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
				m_node_r_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
				m_node_v_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
				m_node_r_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
			});

			res_norm = lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z) / res_norm_0;

			const scalar rho_criterion = (m_pcg_criterion * res_norm_0) * (m_pcg_criterion * res_norm_0);

			scalar rho_old, beta;
			for (; iter < m_maxiters && res_norm > m_pcg_criterion && rho > rho_criterion; ++iter)
			{
				rho_old = rho;

				performInvLocalSolve(scene, m_node_r_x, m_node_r_y, m_node_r_z,
				                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
				                     m_node_z_x, m_node_z_y, m_node_z_z);

				rho = dotNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z, m_node_z_x, m_node_z_y, m_node_z_z);

				beta = rho / rho_old;

				buckets.for_each_bucket([&] (int bucket_idx) {
					m_node_p_x[bucket_idx] = m_node_z_x[bucket_idx] + m_node_p_x[bucket_idx] * beta;
					m_node_p_y[bucket_idx] = m_node_z_y[bucket_idx] + m_node_p_y[bucket_idx] * beta;
					m_node_p_z[bucket_idx] = m_node_z_z[bucket_idx] + m_node_p_z[bucket_idx] * beta;
				});

				performGlobalMultiply(scene, dt,
				                      m_node_Cs_x, m_node_Cs_y, m_node_Cs_z,
				                      m_node_p_x, m_node_p_y, m_node_p_z,
				                      m_node_q_x, m_node_q_y, m_node_q_z);

				alpha = rho / dotNodeVectors(m_node_p_x, m_node_p_y, m_node_p_z, m_node_q_x, m_node_q_y, m_node_q_z);

				buckets.for_each_bucket([&] (int bucket_idx) {
					m_node_v_plus_x[bucket_idx] += m_node_p_x[bucket_idx] * alpha;
					m_node_r_x[bucket_idx] -= m_node_q_x[bucket_idx] * alpha;
					m_node_v_plus_y[bucket_idx] += m_node_p_y[bucket_idx] * alpha;
					m_node_r_y[bucket_idx] -= m_node_q_y[bucket_idx] * alpha;
					m_node_v_plus_z[bucket_idx] += m_node_p_z[bucket_idx] * alpha;
					m_node_r_z[bucket_idx] -= m_node_q_z[bucket_idx] * alpha;
				});

				res_norm = lengthNodeVectors(m_node_r_x, m_node_r_y, m_node_r_z) / res_norm_0;

				if (scene.getLiquidInfo().iteration_print_step > 0 && iter % scene.getLiquidInfo().iteration_print_step == 0)
					std::cout << "[pcg total iter: " << iter
					          << ", res: " << res_norm << "/" << m_pcg_criterion
					          << ", abs. res: " << (res_norm * res_norm_0) << "/" << (m_pcg_criterion * res_norm_0)
					          << ", rho: " << (rho / (res_norm_0 * res_norm_0)) << "/" << (rho_criterion / (res_norm_0 * res_norm_0))
					          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;

			}

			std::cout << "[pcg total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_0) << "/" << (m_pcg_criterion * res_norm_0)
			          << ", rho: " << (rho / (res_norm_0 * res_norm_0)) << "/" << (rho_criterion / (res_norm_0 * res_norm_0))
			          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;
		}
	}

	if (res_norm_1 > m_pcg_criterion)
	{
		const int num_elasto = scene.getNumSoftElastoParticles();

		constructAngularHessianPreProcess(scene, dt);
		constructAngularHessianPostProcess(scene, dt);

		m_angular_r.resize(num_elasto);
		m_angular_z.resize(num_elasto);
		m_angular_p.resize(num_elasto);
		m_angular_q.resize(num_elasto);

		m_angular_r.setZero();
		m_angular_z.setZero();
		m_angular_p.setZero();
		m_angular_q.setZero();

		performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_v_plus_buffer, m_angular_r);

		m_angular_r = m_angular_moment_buffer - m_angular_r;

		scalar res_norm = m_angular_r.norm() / res_norm_1;

		int iter = 0;

		if (res_norm < m_pcg_criterion) {
			std::cout << "[angular pcg total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_1) << "/" << (m_pcg_criterion * res_norm_1)
			          << "]" << std::endl;
		} else {
			performLocalSolveTwist(scene, m_angular_r, scene.getM(), m_angular_z);

			m_angular_p = m_angular_z;

			performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_p, m_angular_q);

			scalar rho = m_angular_r.dot(m_angular_z);

			scalar alpha = rho / m_angular_p.dot(m_angular_q);

			m_angular_v_plus_buffer += m_angular_p * alpha;
			m_angular_r -= m_angular_q * alpha;

			res_norm = m_angular_r.norm() / res_norm_1;

			const scalar rho_criterion = (m_pcg_criterion * res_norm_1) * (m_pcg_criterion * res_norm_1);

			scalar rho_old, beta;
			for (; iter < m_maxiters && res_norm > m_pcg_criterion && rho > rho_criterion; ++iter)
			{
				rho_old = rho;

				performLocalSolveTwist(scene, m_angular_r, scene.getM(), m_angular_z);

				rho = m_angular_r.dot(m_angular_z);

				beta = rho / rho_old;

				m_angular_p = m_angular_z + m_angular_p * beta;

				performAngularGlobalMultiply(scene, dt, scene.getM(), m_angular_p, m_angular_q);

				alpha = rho / m_angular_p.dot(m_angular_q);

				m_angular_v_plus_buffer += m_angular_p * alpha;
				m_angular_r -= m_angular_q * alpha;

				res_norm = m_angular_r.norm() / res_norm_1;

				if (scene.getLiquidInfo().iteration_print_step > 0 && iter % scene.getLiquidInfo().iteration_print_step == 0)
					std::cout << "[angular pcg total iter: " << iter
					          << ", res: " << res_norm << "/" << m_pcg_criterion
					          << ", abs. res: " << (res_norm * res_norm_1) << "/" << (m_pcg_criterion * res_norm_1)
					          << ", rho: " << (rho / (res_norm_1 * res_norm_1)) << "/" << (rho_criterion / (res_norm_1 * res_norm_1))
					          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;

			}

			std::cout << "[angular pcg total iter: " << iter
			          << ", res: " << res_norm << "/" << m_pcg_criterion
			          << ", abs. res: " << (res_norm * res_norm_1) << "/" << (m_pcg_criterion * res_norm_1)
			          << ", rho: " << (rho / (res_norm_1 * res_norm_1)) << "/" << (rho_criterion / (res_norm_1 * res_norm_1))
			          << ", abs. rho: " << rho << "/" << rho_criterion << "]" << std::endl;
		}
	}
	return true;
}

bool LinearizedImplicitEuler::stepImplicitViscosityDiagonalPCG( const TwoDScene& scene,
        const std::vector< VectorXs >& node_vel_src_x,
        const std::vector< VectorXs >& node_vel_src_y,
        const std::vector< VectorXs >& node_vel_src_z,
        std::vector< VectorXs >& node_vel_x,
        std::vector< VectorXs >& node_vel_y,
        std::vector< VectorXs >& node_vel_z,
        const scalar& dt )
{
	allocateNodeVectors(scene, m_node_visc_indices_x, m_node_visc_indices_y, m_node_visc_indices_z);

	int offset_nodes_x;
	int offset_nodes_y;
	int offset_nodes_z;

	const scalar sub_dt = dt / (scalar) m_viscosity_substeps;

	for (int i = 0; i < m_viscosity_substeps; ++i)
	{
		if (i == 0) {
			viscosity::constructViscosityMatrixRHS(scene, m_node_visc_indices_x, m_node_visc_indices_y, m_node_visc_indices_z,
			                                       m_effective_node_indices_x, m_effective_node_indices_y, m_effective_node_indices_z,
			                                       node_vel_src_x, node_vel_src_y, node_vel_src_z,
			                                       m_visc_matrix, m_visc_rhs,
			                                       offset_nodes_x, offset_nodes_y, offset_nodes_z, sub_dt);
		} else {
			viscosity::updateViscosityRHS(scene, m_node_visc_indices_x, m_node_visc_indices_y, m_node_visc_indices_z,
			                              m_effective_node_indices_x, m_effective_node_indices_y, m_effective_node_indices_z,
			                              node_vel_src_x, node_vel_src_y, node_vel_src_z, m_visc_rhs,
			                              offset_nodes_x, offset_nodes_y, offset_nodes_z, sub_dt);
		}


		int iter_out;
		scalar residual;

		viscosity::applyNodeViscosityImplicit(scene, m_node_visc_indices_x, m_node_visc_indices_y, m_node_visc_indices_z,
		                                      offset_nodes_x, offset_nodes_y, offset_nodes_z,
		                                      m_visc_matrix, m_visc_rhs, m_visc_solution,
		                                      node_vel_x, node_vel_y, node_vel_z,
		                                      residual, iter_out,
		                                      m_viscous_criterion, m_maxiters);

		std::cout << "[implicit viscosity sub-step: " << i << ", total iter: " << iter_out << ", res: " << residual << "]" << std::endl;
	}



	return true;
}

bool LinearizedImplicitEuler::acceptVelocity( TwoDScene& scene )
{
	const Sorter& buckets = scene.getParticleBuckets();

	if (scene.getLiquidInfo().solve_solid && scene.getNumSoftElastoParticles() > 0) {
		buckets.for_each_bucket([&] (int bucket_idx) {
			scene.getNodeVelocityX()[bucket_idx] = m_node_v_plus_x[bucket_idx];
			scene.getNodeVelocityY()[bucket_idx] = m_node_v_plus_y[bucket_idx];
			scene.getNodeVelocityZ()[bucket_idx] = m_node_v_plus_z[bucket_idx];
		});

		const int num_elasto = scene.getNumSoftElastoParticles();
		threadutils::for_each(0, num_elasto, [&] (int pidx) {
			scene.getV()[pidx * 4 + 3] = m_angular_v_plus_buffer[pidx];
		});
	}

	if (scene.getNumFluidParticles() > 0) {
		buckets.for_each_bucket([&] (int bucket_idx) {
			scene.getNodeFluidVelocityX()[bucket_idx] = m_node_v_fluid_plus_x[bucket_idx];
			scene.getNodeFluidVelocityY()[bucket_idx] = m_node_v_fluid_plus_y[bucket_idx];
			scene.getNodeFluidVelocityZ()[bucket_idx] = m_node_v_fluid_plus_z[bucket_idx];
		});
	}

	return true;
}

bool LinearizedImplicitEuler::stepImplicitElastoAMGPCG( TwoDScene& scene, scalar dt )
{
	int ndof_elasto = scene.getNumSoftElastoParticles() * 4;

	if (ndof_elasto == 0) return true;

	scalar res_norm_0 = lengthNodeVectors(m_node_rhs_x, m_node_rhs_y, m_node_rhs_z);

	if (res_norm_0 > m_pcg_criterion) {
		// construct Particle Hessian
		constructHessianPreProcess(scene, dt);
		m_A.resize(ndof_elasto, ndof_elasto);
		m_A.reserve(m_triA.size());
		m_A.setFromTriplets(m_triA.begin(), m_triA.end());

		// construct W matrix
		buildLocalGlobalMapping( scene,
		                         m_node_global_indices_x,
		                         m_node_global_indices_y,
		                         m_node_global_indices_z,
		                         m_effective_node_indices,
		                         m_dof_ijk);

		buildNodeToSoftParticlesMat( scene,
		                             m_node_global_indices_x,
		                             m_node_global_indices_y,
		                             m_node_global_indices_z,
		                             m_effective_node_indices,
		                             m_tri_W,
		                             m_W);

		m_A = m_W.transpose() * m_A * m_W;

		const int system_size = m_effective_node_indices.size();

		// finalize LHS
		m_H.resize(system_size);
		m_H.zero();

		// Here we actually used A^T, but this doesn't matter since A should be symmetric
		threadutils::for_each(0, (int) m_A.outerSize(), [&] (int row_idx) {
			for (SparseXs::InnerIterator it(m_A, row_idx); it; ++it)
			{
				int col_idx = it.index();
				m_H.add_to_element(row_idx, col_idx, it.value() * dt * dt);
			}

			const Vector3i& index = m_effective_node_indices[row_idx];
			switch (index(1)) {
			case 0:
				m_H.add_to_element(row_idx, row_idx, m_node_Cs_x[index[0]][index[2]]);
				break;
			case 1:
				m_H.add_to_element(row_idx, row_idx, m_node_Cs_y[index[0]][index[2]]);
				break;
			case 2:
				m_H.add_to_element(row_idx, row_idx, m_node_Cs_z[index[0]][index[2]]);
				break;
			default:
				m_H.add_to_element(row_idx, row_idx, 1.0);
				break;
			}
		});

		// construct RHS
		m_elasto_rhs.resize(system_size);
		m_elasto_result.resize(system_size);

		threadutils::for_each(0, system_size, [&] (int nidx) {
			const Vector3i& index = m_effective_node_indices[nidx];
			switch (index(1)) {
			case 0:
				m_elasto_rhs[nidx] = m_node_rhs_x[index[0]][index[2]];
				m_elasto_result[nidx] = 0.0;//node_vel_x[index[0]][index[2]];
				break;
			case 1:
				m_elasto_rhs[nidx] = m_node_rhs_y[index[0]][index[2]];
				m_elasto_result[nidx] = 0.0;//node_vel_y[index[0]][index[2]];
				break;
			case 2:
				m_elasto_rhs[nidx] = m_node_rhs_z[index[0]][index[2]];
				m_elasto_result[nidx] = 0.0;//node_vel_z[index[0]][index[2]];
				break;
			default:
				m_elasto_rhs[nidx] = 0.0;
				m_elasto_result[nidx] = 0.0;
				break;
			}
		});

		const Sorter& buckets = scene.getParticleBuckets();
		const int bucket_num_cell = scene.getDefaultNumNodes();
		const int ni = buckets.ni * bucket_num_cell;
		const int nj = buckets.nj * bucket_num_cell;
		const int nk = buckets.nk * bucket_num_cell;

		bool success = false;
		scalar tolerance = 0.0;
		int iterations = 0;

		success = AMGPCGSolveSparse(m_H, m_elasto_rhs, m_elasto_result,
		                            m_dof_ijk, m_pcg_criterion, m_maxiters,
		                            tolerance, iterations, ni * 3, nj, nk);

		std::cout << "[amg pcg elasto total iter: " << iterations << ", res: " << tolerance << "]" << std::endl;

		if (!success) {
			std::cout << "WARNING: AMG PCG solve failed!" << std::endl;

			std::cout << "rhs=[";
			for (scalar s : m_elasto_rhs) {
				std::cout << s << "; ";
			}
			std::cout << "];" << std::endl;
		}

		threadutils::for_each(0, system_size, [&] (int nidx) {
			const Vector3i& index = m_effective_node_indices[nidx];
			switch (index(1)) {
			case 0:
				m_node_v_plus_x[index[0]][index[2]] = m_elasto_result[nidx];
				break;
			case 1:
				m_node_v_plus_y[index[0]][index[2]] = m_elasto_result[nidx];
				break;
			case 2:
				m_node_v_plus_z[index[0]][index[2]] = m_elasto_result[nidx];
				break;
			default:
				break;
			}
		});
	}

	return true;
}

bool LinearizedImplicitEuler::stepImplicitElasto( TwoDScene& scene, scalar dt )
{
	if (scene.getLiquidInfo().use_amgpcg_solid) {
		return stepImplicitElastoAMGPCG(scene, dt);
	} else if (scene.getLiquidInfo().use_pcr) {
		return stepImplicitElastoDiagonalPCR(scene, dt);
	} else {
		if (scene.getLiquidInfo().use_cosolve_angular)
			return stepImplicitElastoDiagonalPCGCoSolve(scene, dt);
		else
			return stepImplicitElastoDiagonalPCG(scene, dt);
	}
}

bool LinearizedImplicitEuler::applyPressureDragElasto( TwoDScene& scene, scalar dt )
{
	if (scene.getNumFluidParticles() == 0) return false;

	int ndof_elasto = scene.getNumSoftElastoParticles() * 4;

	// u_f^*
	const auto& m_node_v_f_star_x = scene.getNodeFluidVelocityX();
	const auto& m_node_v_f_star_y = scene.getNodeFluidVelocityY();
	const auto& m_node_v_f_star_z = scene.getNodeFluidVelocityZ();

	// solve elasto
	if (ndof_elasto > 0) {
		if (scene.getLiquidInfo().apply_pressure_solid)
			pressure::applyPressureGradsElastoRHS(scene, scene.getNodePressure(), m_node_rhs_x, m_node_rhs_y, m_node_rhs_z, m_node_mfhdvm_hdvm_x, m_node_mfhdvm_hdvm_y, m_node_mfhdvm_hdvm_z, dt);
		//
		addFluidDragRHS(scene, m_node_v_f_star_x, m_node_v_f_star_y, m_node_v_f_star_z, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z);

		performInvLocalSolve(scene, m_node_rhs_x, m_node_rhs_y, m_node_rhs_z,
		                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
		                     m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z);
	}

	return true;
}

void LinearizedImplicitEuler::pushFluidVelocity()
{
	m_fluid_vel_stack.push(m_node_v_fluid_plus_x);
	m_fluid_vel_stack.push(m_node_v_fluid_plus_y);
	m_fluid_vel_stack.push(m_node_v_fluid_plus_z);

	m_fluid_vel_stack.push(m_node_rhs_fluid_x);
	m_fluid_vel_stack.push(m_node_rhs_fluid_y);
	m_fluid_vel_stack.push(m_node_rhs_fluid_z);
}

void LinearizedImplicitEuler::popFluidVelocity()
{
	m_node_rhs_fluid_z = m_fluid_vel_stack.top();
	m_fluid_vel_stack.pop();
	m_node_rhs_fluid_y = m_fluid_vel_stack.top();
	m_fluid_vel_stack.pop();
	m_node_rhs_fluid_x = m_fluid_vel_stack.top();
	m_fluid_vel_stack.pop();

	m_node_v_fluid_plus_z = m_fluid_vel_stack.top();
	m_fluid_vel_stack.pop();
	m_node_v_fluid_plus_y = m_fluid_vel_stack.top();
	m_fluid_vel_stack.pop();
	m_node_v_fluid_plus_x = m_fluid_vel_stack.top();
	m_fluid_vel_stack.pop();
}

void LinearizedImplicitEuler::pushElastoVelocity()
{
	m_elasto_vel_stack.push(m_node_v_plus_x);
	m_elasto_vel_stack.push(m_node_v_plus_y);
	m_elasto_vel_stack.push(m_node_v_plus_z);

	m_elasto_vel_stack.push(m_node_rhs_x);
	m_elasto_vel_stack.push(m_node_rhs_y);
	m_elasto_vel_stack.push(m_node_rhs_z);
}

void LinearizedImplicitEuler::popElastoVelocity()
{
	m_node_rhs_z = m_elasto_vel_stack.top();
	m_elasto_vel_stack.pop();
	m_node_rhs_y = m_elasto_vel_stack.top();
	m_elasto_vel_stack.pop();
	m_node_rhs_x = m_elasto_vel_stack.top();
	m_elasto_vel_stack.pop();

	m_node_v_plus_z = m_elasto_vel_stack.top();
	m_elasto_vel_stack.pop();
	m_node_v_plus_y = m_elasto_vel_stack.top();
	m_elasto_vel_stack.pop();
	m_node_v_plus_x = m_elasto_vel_stack.top();
	m_elasto_vel_stack.pop();
}

bool LinearizedImplicitEuler::applyPressureDragFluid( TwoDScene& scene, scalar dt )
{
	if (scene.getNumFluidParticles() == 0) return false;

	const std::vector< VectorXs >& node_mass_fluid_x = scene.getNodeFluidMassX();
	const std::vector< VectorXs >& node_mass_fluid_y = scene.getNodeFluidMassY();
	const std::vector< VectorXs >& node_mass_fluid_z = scene.getNodeFluidMassZ();

	// update liquid velocity
	if (scene.getNumFluidParticles() > 0) {
		// apply pressure
		std::vector< VectorXs >& node_pressure = scene.getNodePressure();

		if (scene.getLiquidInfo().drag_by_future_solid) {
			const Sorter& buckets = scene.getParticleBuckets();

			buckets.for_each_bucket([&] (int bucket_idx) {
				const int num_nodes = scene.getNumNodes(bucket_idx);

				for (int i = 0; i < num_nodes; ++i)
				{
					m_node_v_fluid_plus_x[bucket_idx][i] *= node_mass_fluid_x[bucket_idx][i] * m_node_inv_mfhdvm_x[bucket_idx][i];
					m_node_v_fluid_plus_y[bucket_idx][i] *= node_mass_fluid_y[bucket_idx][i] * m_node_inv_mfhdvm_y[bucket_idx][i];
					m_node_v_fluid_plus_z[bucket_idx][i] *= node_mass_fluid_z[bucket_idx][i] * m_node_inv_mfhdvm_z[bucket_idx][i];
				}
			});

			addSolidDrag(scene, m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z);

			pressure::applyPressureGradsFluid(scene, node_pressure, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z, m_node_inv_mfhdvm_x, m_node_inv_mfhdvm_y, m_node_inv_mfhdvm_z, dt);
		} else {
			const auto& m_node_u_s_star_x = scene.getNodeVelocityX();
			const auto& m_node_u_s_star_y = scene.getNodeVelocityY();
			const auto& m_node_u_s_star_z = scene.getNodeVelocityZ();

			addSolidDragRHS(scene, m_node_u_s_star_x, m_node_u_s_star_y, m_node_u_s_star_z, m_node_rhs_fluid_x, m_node_rhs_fluid_y, m_node_rhs_fluid_z);

			pressure::applyPressureGradsFluidRHS(scene, scene.getNodePressure(), m_node_rhs_fluid_x, m_node_rhs_fluid_y, m_node_rhs_fluid_z, m_node_mshdvm_hdvm_x, m_node_mshdvm_hdvm_y, m_node_mshdvm_hdvm_z, dt);

			performInvLocalSolve(scene, m_node_rhs_fluid_x, m_node_rhs_fluid_y, m_node_rhs_fluid_z, m_node_inv_C_x, m_node_inv_C_y, m_node_inv_C_z, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z);
		}

		pressure::extrapolate(scene, scene.getNodeLiquidValidX(), m_node_v_fluid_plus_x);
		pressure::extrapolate(scene, scene.getNodeLiquidValidY(), m_node_v_fluid_plus_y);
		pressure::extrapolate(scene, scene.getNodeLiquidValidZ(), m_node_v_fluid_plus_z);
	}

	return true;
}

bool LinearizedImplicitEuler::projectFine( TwoDScene& scene, scalar dt )
{
	if (scene.getNumFluidParticles() == 0) return false;

	allocateCenterNodeVectors(scene, m_fine_global_indices);

	pressure::solveNodePressure(scene, scene.getNodePressure(), m_fine_pressure_rhs,
	                            m_fine_pressure_matrix, m_fine_global_indices,
	                            m_node_psi_fs_x, m_node_psi_fs_y, m_node_psi_fs_z,
	                            m_node_psi_sf_x, m_node_psi_sf_y, m_node_psi_sf_z,
	                            m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z,
	                            m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z,
	                            m_node_inv_C_x, m_node_inv_C_y, m_node_inv_C_z,
	                            m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
	                            m_node_mfhdvm_hdvm_x, m_node_mfhdvm_hdvm_y, m_node_mfhdvm_hdvm_z,
	                            m_node_mshdvm_hdvm_x, m_node_mshdvm_hdvm_y, m_node_mshdvm_hdvm_z,
	                            dt, m_pressure_criterion, m_maxiters);

#ifdef CHECK_EQU_24
	pushFluidVelocity();
	pushElastoVelocity();

	// apply pressure & check divergence on the Naive Equ.
	// construct u_s^n+1 according to Equ. 24
	std::vector< VectorXs > tmp_rhs_x = m_node_rhs_x;
	std::vector< VectorXs > tmp_rhs_y = m_node_rhs_y;
	std::vector< VectorXs > tmp_rhs_z = m_node_rhs_z;

	const auto& m_node_u_f_star_x = scene.getNodeFluidVelocityX();
	const auto& m_node_u_f_star_y = scene.getNodeFluidVelocityY();
	const auto& m_node_u_f_star_z = scene.getNodeFluidVelocityZ();

	int num_elasto = scene.getNumSoftElastoParticles();
	// solve elasto
	if (num_elasto > 0) {
		if (scene.getLiquidInfo().apply_pressure_solid)
			pressure::applyPressureGradsElastoRHS(scene, scene.getNodePressure(), tmp_rhs_x, tmp_rhs_y, tmp_rhs_z, m_node_mfhdvm_hdvm_x, m_node_mfhdvm_hdvm_y, m_node_mfhdvm_hdvm_z, dt);

		addFluidDragRHS(scene, m_node_u_f_star_x, m_node_u_f_star_y, m_node_u_f_star_z, tmp_rhs_x, tmp_rhs_y, tmp_rhs_z);

		performInvLocalSolve(scene, tmp_rhs_x, tmp_rhs_y, tmp_rhs_z,
		                     m_node_inv_Cs_x, m_node_inv_Cs_y, m_node_inv_Cs_z,
		                     m_node_v_plus_x, m_node_v_plus_y, m_node_v_plus_z);
	}

	// construct u_f^n+1 according to Equ. 24
	std::vector< VectorXs > tmp_rhs_fluid_x = m_node_rhs_fluid_x;
	std::vector< VectorXs > tmp_rhs_fluid_y = m_node_rhs_fluid_y;
	std::vector< VectorXs > tmp_rhs_fluid_z = m_node_rhs_fluid_z;

	const auto& m_node_u_s_star_x = scene.getNodeVelocityX();
	const auto& m_node_u_s_star_y = scene.getNodeVelocityY();
	const auto& m_node_u_s_star_z = scene.getNodeVelocityZ();

	int num_fluids = scene.getNumFluidParticles();
	// solve fluid
	if (num_fluids > 0) {
		addSolidDragRHS(scene, m_node_u_s_star_x, m_node_u_s_star_y, m_node_u_s_star_z, tmp_rhs_fluid_x, tmp_rhs_fluid_y, tmp_rhs_fluid_z);

		pressure::applyPressureGradsFluidRHS(scene, scene.getNodePressure(), tmp_rhs_fluid_x, tmp_rhs_fluid_y, tmp_rhs_fluid_z, m_node_mshdvm_hdvm_x, m_node_mshdvm_hdvm_y, m_node_mshdvm_hdvm_z, dt);

		performInvLocalSolve(scene, tmp_rhs_fluid_x, tmp_rhs_fluid_y, tmp_rhs_fluid_z, m_node_inv_C_x, m_node_inv_C_y, m_node_inv_C_z, m_node_v_fluid_plus_x, m_node_v_fluid_plus_y, m_node_v_fluid_plus_z);

		pressure::extrapolate(scene, scene.getNodeLiquidValidX(), scene.getNodeCompressedIndexX(), m_node_v_fluid_plus_x);
		pressure::extrapolate(scene, scene.getNodeLiquidValidY(), scene.getNodeCompressedIndexY(), m_node_v_fluid_plus_y);
		pressure::extrapolate(scene, scene.getNodeLiquidValidZ(), scene.getNodeCompressedIndexZ(), m_node_v_fluid_plus_z);
	}

	// check divergence
	scalar div = computeDivergence(scene);

	std::cout << "CHECK EQU 24: " << div << std::endl;

	popFluidVelocity();
	popElastoVelocity();
#endif
	return true;
}

std::string LinearizedImplicitEuler::getName() const
{
	return "Linearized Implicit Euler";
}

void LinearizedImplicitEuler::zeroFixedDoFs( const TwoDScene& scene, VectorXs& vec )
{
	int nprts = scene.getNumParticles();
	threadutils::for_each(0, nprts, [&] (int i) {
		if ( scene.isFixed(i) & 1 ) vec.segment<3>(4 * i).setZero();
		if ( scene.isFixed(i) & 2 ) vec(4 * i + 3) = 0.0;
	});
}

bool LinearizedImplicitEuler::manifoldPropagate( TwoDScene& scene, scalar dt )
{
	const scalar subdt = dt / (scalar) m_manifold_substeps;

	VectorXs& fluid_vol = scene.getFluidVol();
	VectorXs& fluid_m = scene.getFluidM();
	VectorXs& elasto_v = scene.getV();
	VectorXs& elasto_m = scene.getM();

	const VectorXs& fluid_m_gauss = scene.getGaussFluidM();
	const VectorXs& fluid_vol_gauss = scene.getGaussFluidVol();
	const VectorXs& vol_gauss = scene.getGaussVol();
	const VectorXs& vol_frac_gauss = scene.getGaussVolumeFraction();

	const VectorXs& dv_gauss = scene.getGaussDV();
	const VectorXs& v_gauss = scene.getGaussV();
	const std::vector< VectorXs >& part_div = scene.getParticleDiv();
	const std::vector< int > particle_to_surfels = scene.getParticleToSurfels();

	int total_iter = 0;
	scalar total_res = 0.0;
	for (int k = 0; k < m_manifold_substeps; ++k) {
		const int num_gauss = scene.getNumGausses();
		VectorXs F(num_gauss * 3);

		const int num_edges = scene.getNumEdges();
		const int num_faces = scene.getNumFaces();

		const int num_elasto_gauss = num_edges + num_faces;
		const int num_elasto = scene.getNumSoftElastoParticles();

		VectorXs fv0 = fluid_vol;

		scalar res_0 = fv0.segment(0, num_elasto).norm();

		scalar res_norm = 1.0;

		if (res_0 < 1e-20) res_norm = 0.0;


		int iter = 0;
		for (; iter < m_maxiters && res_norm > m_quasi_static_criterion; ++iter)
		{
			VectorXs old_fv = fluid_vol;

			F.setZero();
			scene.accumulateManifoldFluidGradU(F);
			scene.accumulateManifoldGradPorePressure(F);

			threadutils::for_each(0, num_elasto_gauss, [&] (int gidx) {
				F.segment<3>(gidx * 3) -= dv_gauss.segment<3>(gidx * 4) / dt * fluid_m_gauss(gidx * 4);

				F.segment<3>(gidx * 3) *= (1.0 - vol_frac_gauss[gidx]);

				Vector3s dv = Vector3s::Zero();
				if (fluid_m_gauss(gidx * 4) > 1e-20)
					dv = F.segment<3>(gidx * 3) / fluid_m_gauss(gidx * 4) * dt;

				const scalar min_D = fluid_m_gauss(gidx * 4) / dt * 1e-3;
				const scalar vol_empty = vol_gauss(gidx) * (1.0 - vol_frac_gauss(gidx));
				const scalar s = (vol_empty > 1e-20) ? mathutils::clamp(fluid_vol_gauss(gidx) / vol_empty, 0.0, 1.0) : 0.0;

				const scalar D = std::max(min_D, scene.getPlanarDragCoeff(vol_frac_gauss(gidx), s, dv.norm(), 0));

				if (D > 1e-20) F.segment<3>(gidx * 3) /= D;
			});

			if (scene.propagateSolidVelocity()) {
				threadutils::for_each(0, num_elasto, [&] (int pidx) {
					if (particle_to_surfels[pidx] >= 0) return;

					const VectorXs& div = part_div[pidx];
					const auto& pe = scene.getParticleEdges(pidx);
					const int num_pe = pe.size();

					Vector3s divFV = Vector3s::Zero();
					scalar divF = 0.0;
					scalar w = 0.0;

					for (int j = 0; j < num_pe; ++j) {
						const int eidx = pe[j];

						const scalar coeff = F.segment<3>(eidx * 3).dot(div.segment<3>(j * 3)) * fluid_vol_gauss(eidx);
						divF += coeff;
						divFV += v_gauss.segment<3>(eidx * 4) * coeff;
						w += fluid_vol_gauss(eidx);
					}

					const auto& pf = scene.getParticleFaces(pidx);
					const int num_pf = pf.size();

					for (int j = 0; j < num_pf; ++j) {
						const int fidx = pf[j].first;
						const int gidx = fidx + num_edges;

						const scalar coeff = F.segment<3>(gidx * 3).dot(div.segment<3>((j + num_pe) * 3)) * fluid_vol_gauss(gidx);
						divF += coeff;
						divFV += v_gauss.segment<3>(gidx * 4) * coeff;
						w += fluid_vol_gauss(gidx);
					}

					if (w > 1e-20) {
						divF /= w;
						divFV /= w;
					}

					const scalar old_m = elasto_m(pidx * 4) + fluid_m(pidx * 4);

					const scalar new_fluid_vol = std::max(0.0, fv0(pidx) - subdt * divF);

					fluid_vol(pidx) = new_fluid_vol;

					const scalar new_fluid_m = new_fluid_vol * scene.getLiquidInfo().liquid_density;

					const scalar new_m = elasto_m(pidx * 4) + new_fluid_m;

					fluid_m.segment<3>(pidx * 4).setConstant( new_fluid_m );

					const Vector3s old_inertia = old_m * elasto_v.segment<3>(pidx * 4);

					const Vector3s new_inertia = old_inertia - subdt * divFV * scene.getLiquidInfo().liquid_density;

					if (new_m > 1e-20) {
						elasto_v.segment<3>(pidx * 4) = new_inertia / new_m;
					}
				});
			} else {
				threadutils::for_each(0, num_elasto, [&] (int pidx) {
					if (particle_to_surfels[pidx] >= 0) return;

					const VectorXs& div = part_div[pidx];
					const auto& pe = scene.getParticleEdges(pidx);
					const int num_pe = pe.size();

					scalar divF = 0.0;
					scalar w = 0.0;

					for (int j = 0; j < num_pe; ++j) {
						const int eidx = pe[j];

						const scalar coeff = F.segment<3>(eidx * 3).dot(div.segment<3>(j * 3)) * fluid_vol_gauss(eidx);
						divF += coeff;
						w += fluid_vol_gauss(eidx);
					}

					const auto& pf = scene.getParticleFaces(pidx);
					const int num_pf = pf.size();

					for (int j = 0; j < num_pf; ++j) {
						const int fidx = pf[j].first;
						const int gidx = fidx + num_edges;

						const scalar coeff = F.segment<3>(gidx * 3).dot(div.segment<3>((j + num_pe) * 3)) * fluid_vol_gauss(gidx);
						divF += coeff;
						w += fluid_vol_gauss(gidx);
					}

					if (w > 1e-20) {
						divF /= w;
					}

					const scalar old_m = elasto_m(pidx * 4) + fluid_m(pidx * 4);

					const scalar new_fluid_vol = std::max(0.0, fv0(pidx) - subdt * divF);

					fluid_vol(pidx) = new_fluid_vol;

					const scalar new_fluid_m = new_fluid_vol * scene.getLiquidInfo().liquid_density;

					const scalar new_m = elasto_m(pidx * 4) + new_fluid_m;

					fluid_m.segment<3>(pidx * 4).setConstant( new_fluid_m );

					if (new_m > 1e-20) {
						const scalar prop = mathutils::clamp(old_m / new_m, 0.0, 1.0);
						elasto_v.segment<4>(pidx * 4) *= prop;
					}
				});
			}

			scalar old_sum_fv = old_fv.segment(0, num_elasto).sum();
			if (old_sum_fv > 1e-20) {
				scalar new_sum_fv = fluid_vol.segment(0, num_elasto).sum();
				scalar prop = std::min(1.0, old_sum_fv / new_sum_fv);
				fluid_vol.segment(0, num_elasto) *= prop;
				fluid_m.segment(0, num_elasto * 4) *= prop;
			}


			scene.updateGaussManifoldSystem();
			scene.updateVelocityDifference();
			scene.updateGaussAccel();

			res_norm = (fluid_vol.segment(0, num_elasto) - old_fv.segment(0, num_elasto)).norm() / res_0;
		}

		total_iter += iter;
		total_res += res_norm * res_norm;
	}

	std::cout << "[manifold propagate avg iter: " << ((scalar) total_iter / (scalar) m_manifold_substeps) << ", avg res: " << sqrt(total_res / (scalar) m_manifold_substeps) << "]" << std::endl;

	return true;
}
