//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "JunctionForce.h"
#include "TwoDScene.h"

JunctionForce::JunctionForce(const std::shared_ptr< TwoDScene >& scene)
	: m_scene(scene)
{
	const int num_elasto = scene->getNumSoftElastoParticles();
	const int num_edges = scene->getNumEdges();
	const VectorXs& radius = scene->getRadius();
	const MatrixXi& edges = scene->getEdges();

	m_junctions_indices.reserve(num_elasto);
	m_junctions_edges.reserve(num_elasto);
	m_base_indices.reserve(num_elasto);
	m_bending_coeff.reserve(num_elasto);

	m_count_edges = 0;

	for (int pidx = 0; pidx < num_elasto; ++pidx) {
		if (scene->getParticleEdges(pidx).size() > 0 && scene->getParticleFaces(pidx).size() > 0) {
			m_junctions_indices.push_back(pidx);

			const auto& pes = scene->getParticleEdges(pidx);

			std::vector<int> another_ps;
			std::vector<scalar> another_coeff;
			for (int e : pes) {
				int another_p = (edges(e, 0) == pidx) ? edges(e, 1) : edges(e, 0);
				another_ps.push_back(another_p);

				const scalar ra = radius(another_p * 2 + 0);
				const scalar rb = radius(another_p * 2 + 1);
				const scalar coeff = M_PI / 4.0 * (ra * rb * (ra * ra + rb * rb)) * m_scene->getYoungModulus(e) / m_scene->getParticleRestLength(another_p);
				another_coeff.push_back(coeff);
			}

			m_bending_coeff.push_back(another_coeff);
			m_junctions_edges.push_back(another_ps);

			m_base_indices.push_back(m_count_edges);
			m_count_edges += (int) another_ps.size();
		}
	};

	m_junction_orientation.resize(m_junctions_indices.size() * 3);
	m_junction_orientation.setZero();

	preCompute();

	const VectorXs& x = scene->getRestPos();

	m_junction_signs.resize(m_junctions_indices.size());

	threadutils::for_each(0, (int) m_junctions_indices.size(), [&] (int i) {
		const int pidx = m_junctions_indices[i];

		const std::vector<int>& another_ps = m_junctions_edges[i];
		m_junction_signs[i].resize(another_ps.size());

		for (int j = 0; j < (int) another_ps.size(); ++j) {
			Vector3s e = x.segment<3>(another_ps[j] * 4) - x.segment<3>(pidx * 4);
			Vector3s dir = m_junction_orientation.segment<3>(i * 3);

			m_junction_signs[i][j] = mathutils::sgn( e.dot(dir) ) * e.norm();
		}
	});
}

void JunctionForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, scalar& E )
{
	std::cerr << "NOT IMPLEMENTED!" << std::endl;
}

void JunctionForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, VectorXs& gradE )
{
	threadutils::for_each(0, (int) m_junctions_indices.size(), [&] (int i) {
		const int pidx = m_junctions_indices[i];
		const scalar psi_coeff = pow(psi(pidx), lambda);
		const std::vector<int>& another_ps = m_junctions_edges[i];
		const int num_pes = (int) another_ps.size();

		const Vector3s& x0 = x.segment<3>(pidx * 4);
		const Vector3s& ori = m_junction_orientation.segment<3>(i * 3);

		for (int j = 0; j < num_pes; ++j)
		{
			const int npidx = another_ps[j];
			const Vector3s& x1 = x.segment<3>(npidx * 4);

			const Vector3s e12 = x1 - x0;
			const Vector3s e23 = ori * m_junction_signs[i][j];

			const scalar e12n = e12.norm();
			const scalar e23n = e23.norm();

			scalar ece = e12.cross(e23).norm();
			scalar ede = e12.dot(e23);
			scalar ece2pede2 = ece * ece + ede * ede;
			if (ece2pede2 > 1e-63) {
				Vector6s gradece;
				gradece.segment<3>(0) = -e23.cross(e12.cross(e23));
				gradece.segment<3>(3) = -gradece.segment<3>(0);

				Vector6s gradede;
				gradede.segment<3>(0) = -e23;
				gradede.segment<3>(3) = e23;

				Vector6s gradtheta = (gradece * ede - gradede * ece) / ece2pede2;
				const scalar theta = atan2(ece, ede);

				Vector6s localGradE = m_bending_coeff[i][j] * psi_coeff * theta * gradtheta;

				gradE.segment<3>(3 * pidx) += localGradE.segment<3>(0);
				gradE.segment<3>(3 * npidx) += localGradE.segment<3>(3);
			}
		}
	});
}

bool JunctionForce::parallelized() const
{
	return true;
}

void JunctionForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, const VectorXs& psi, const scalar& lambda, TripletXs& hessE, int hessE_index, const scalar& dt )
{
	threadutils::for_each(0, (int) m_junctions_indices.size(), [&] (int i) {
		const int pidx = m_junctions_indices[i];
		const scalar psi_coeff = pow(psi(pidx), lambda);
		const std::vector<int>& another_ps = m_junctions_edges[i];
		const int num_pes = (int) another_ps.size();

		const Vector3s& x0 = x.segment<3>(pidx * 4);
		const Vector3s& ori = m_junction_orientation.segment<3>(i * 3);


		for (int j = 0; j < num_pes; ++j)
		{
			const int base_idx = (m_base_indices[i] + j) * 6 * 6;

			const int npidx = another_ps[j];
			const Vector3s& x1 = x.segment<3>(npidx * 4);

			const Vector3s e12 = x1 - x0;
			const Vector3s e23 = ori * m_junction_signs[i][j];

			const scalar e12n = e12.norm();
			const scalar e23n = e23.norm();

			scalar ece = e12.cross(e23).norm();
			scalar ede = e12.dot(e23);
			scalar ece2pede2 = ece * ece + ede * ede;
			if (ece2pede2 > 1e-63) {
				Vector6s gradece;
				gradece.segment<3>(0) = -e23.cross(e12.cross(e23));
				gradece.segment<3>(3) = -gradece.segment<3>(0);

				Vector6s gradede;
				gradede.segment<3>(0) = -e23;
				gradede.segment<3>(3) = e23;

				Vector6s gradtheta = (gradece * ede - gradede * ece) / ece2pede2;
				const scalar theta = atan2(ece, ede);

				Matrix3s hessece_comp;
				hessece_comp <<
				             e23(1) * e23(1) + e23(2) * e23(2), -e23(0) * e23(1), -e23(0) * e23(2),
				                 -e23(0) * e23(1), e23(0) * e23(0) + e23(2) * e23(2), -e23(1) * e23(2),
				                 -e23(0) * e23(2), -e23(1) * e23(2), e23(1) * e23(1) + e23(2) * e23(2);

				Matrix6s hessece;
				hessece.block<3, 3>(0, 0) = hessece.block<3, 3>(3, 3) = hessece_comp;
				hessece.block<3, 3>(3, 0) = hessece.block<3, 3>(0, 3) = -hessece_comp;

				Matrix6s hesstheta = ((hessece * ede + gradece * gradede.transpose() - gradede * gradece.transpose()) * ece2pede2 - 2.0 * (gradece * ede - gradede * ece) * (ece * gradece + ede * gradede).transpose()) / (ece2pede2 * ece2pede2);

				Matrix6s localhessE = m_bending_coeff[i][j] * psi_coeff * (gradtheta * gradtheta.transpose() + theta * hesstheta);

				for (int r = 0; r < 3; ++r) for (int s = 0; s < 3; ++s) {
						hessE[hessE_index + base_idx + r * 3 + s] = Triplets(pidx * 4 + r, pidx * 4 + s, localhessE(r, s));
						hessE[hessE_index + base_idx + 9 + r * 3 + s] = Triplets(pidx * 4 + r, npidx * 4 + s, localhessE(r, 3 + s));
						hessE[hessE_index + base_idx + 18 + r * 3 + s] = Triplets(npidx * 4 + r, pidx * 4 + s, localhessE(3 + r, s));
						hessE[hessE_index + base_idx + 27 + r * 3 + s] = Triplets(npidx * 4 + r, npidx * 4 + s, localhessE(3 + r, 3 + s));
					}
			} else {
				for (int r = 0; r < 36; ++r) {
					hessE[hessE_index + base_idx + r] = Triplets(0, 0, 0.0);
				}
			}
		}
	});
}

void JunctionForce::updateMultipliers( const VectorXs& x, const VectorXs& vplus, const VectorXs& m, const VectorXs& psi, const scalar& lambda, const scalar& dt )
{

}

int JunctionForce::numHessX()
{
	return m_count_edges * 6 * 6;
}

void JunctionForce::preCompute()
{
	m_junction_orientation.setZero();
	const int num_junc = m_junctions_indices.size();
	const MatrixXs& norm_gauss = m_scene->getGaussNormal();
	const int num_edges = m_scene->getNumEdges();

	threadutils::for_each(0, num_junc, [&] (int i) {
		const int pidx = m_junctions_indices[i];
		const auto& pairs = m_scene->getParticleFaces(pidx);

		Vector3s n = Vector3s::Zero();
		for (auto& p : pairs)
		{
			const int gidx = p.first + num_edges;
			n += norm_gauss.block<3, 1>(gidx * 3, 2);
		}
		m_junction_orientation.segment<3>(i * 3) = n.normalized();
	});
}

void JunctionForce::updateStartState()
{

}

Force* JunctionForce::createNewCopy()
{
	return new JunctionForce(*this);
}

int JunctionForce::flag() const
{
	return 1;
}
