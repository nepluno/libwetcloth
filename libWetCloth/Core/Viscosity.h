//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
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

#ifndef __VISCOSITY_H__
#define __VISCOSITY_H__

#include <Eigen/Core>

#include "MathDefs.h"
#include "array3.h"
#include "pcgsolver/sparse_matrix.h"

class TwoDScene;

using namespace robertbridson;

namespace viscosity {
	
	void constructViscosityMatrixRHS( const TwoDScene& scene,
									 std::vector< VectorXi >& node_global_indices_x,
									 std::vector< VectorXi >& node_global_indices_y,
									 std::vector< VectorXi >& node_global_indices_z,
									 std::vector< Vector2i >& effective_node_indices_x,
									 std::vector< Vector2i >& effective_node_indices_y,
									 std::vector< Vector2i >& effective_node_indices_z,
									 const std::vector< VectorXs >& node_vel_src_x,
									 const std::vector< VectorXs >& node_vel_src_y,
									 const std::vector< VectorXs >& node_vel_src_z,
									 SparseMatrix< scalar >& matrix,
									 std::vector< scalar >& rhs,
									 int& offset_nodes_x,
									 int& offset_nodes_y,
									 int& offset_nodes_z,
									 const scalar& dt );
	
	void updateViscosityRHS( const TwoDScene& scene,
							const std::vector< VectorXi >& node_global_indices_x,
							const std::vector< VectorXi >& node_global_indices_y,
							const std::vector< VectorXi >& node_global_indices_z,
							const std::vector< Vector2i >& effective_node_indices_x,
							const std::vector< Vector2i >& effective_node_indices_y,
							const std::vector< Vector2i >& effective_node_indices_z,
							const std::vector< VectorXs >& node_vel_src_x,
							const std::vector< VectorXs >& node_vel_src_y,
							const std::vector< VectorXs >& node_vel_src_z,
							std::vector< scalar >& rhs,
							int offset_nodes_x,
							int offset_nodes_y,
							int offset_nodes_z,
							const scalar& dt  );
	
	void applyNodeViscosityImplicit( const TwoDScene& scene,
									const std::vector< VectorXi >& node_global_indices_x,
									const std::vector< VectorXi >& node_global_indices_y,
									const std::vector< VectorXi >& node_global_indices_z,
									int offset_nodes_x,
									int offset_nodes_y,
									int offset_nodes_z,
									const SparseMatrix< scalar >& matrix,
									const std::vector< scalar >& rhs,
									std::vector< scalar >& soln,
									std::vector< VectorXs >& node_vel_x,
									std::vector< VectorXs >& node_vel_y,
									std::vector< VectorXs >& node_vel_z,
									scalar& residual,
									int& iter_out,
									const scalar& criterion,
									int maxiters);
	
	void applyNodeViscosityExplicit( const TwoDScene& scene,
									const std::vector< VectorXs >& node_vel_src_x,
									const std::vector< VectorXs >& node_vel_src_y,
									const std::vector< VectorXs >& node_vel_src_z,
									std::vector< VectorXs >& node_vel_x,
									std::vector< VectorXs >& node_vel_y,
									std::vector< VectorXs >& node_vel_z,
									const scalar& dt );
};

#endif
