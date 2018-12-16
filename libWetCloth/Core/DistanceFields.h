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

#ifndef __DISTANCE_FIELDS_H__
#define __DISTANCE_FIELDS_H__

#include <Eigen/Core>
#include <memory>

#include "MathDefs.h"
#include "SolidMesh.h"
#include "array3.h"
#include "array3_utils.h"

class TwoDScene;

struct DF_SOURCE_DURATION
{
	scalar start;
	scalar end;
    scalar maxvol;
	Vector3s vel;
};

enum DISTANCE_FIELD_USAGE
{
	DFU_SOLID,
	DFU_SOURCE,
    DFU_TERMINATOR,
	
	DFU_COUNT
};

enum DISTANCE_FIELD_TYPE
{
	DFT_SPHERE,
	DFT_BOX,
	DFT_CAPSULE,
    DFT_CYLINDER,
    DFT_FILE,
	
	DFT_UNION,
	DFT_INTERSECT,
	
	DFT_COUNT
};

struct DistanceField
{
	DistanceField(DISTANCE_FIELD_TYPE type_, DISTANCE_FIELD_USAGE usage_, int group_, int params_index_, bool sampled_);
	
	virtual void advance(const scalar& dt) = 0;
	virtual scalar compute_phi(const Vector3s& pos) const = 0;
	virtual scalar compute_phi_vel(const Vector3s& pos, Vector3s& vel) const = 0;
	virtual void sample(const scalar& dx, VectorXs& result, VectorXs& normals);
	virtual bool check_durations(const scalar& cur_time, const scalar& cur_vol, Vector3s& shooting_vel) = 0;
    virtual void resample_mesh(const scalar& dx, VectorXs& result, VectorXs& normals) = 0;
    virtual void resample_internal( const std::shared_ptr< TwoDScene >& parent, const scalar& dx, const VectorXs& exist, VectorXs& additional);
	virtual bool local_bounding_box(Vector3s& low, Vector3s& high) const = 0;
	virtual void apply_global_rotation(const Eigen::Quaternion<scalar>& rot) = 0;
	virtual void apply_local_rotation(const Eigen::Quaternion<scalar>& rot) = 0;
	virtual void apply_translation(const Vector3s& t) = 0;
	virtual int vote_param_indices() { return params_index; };
	virtual DISTANCE_FIELD_USAGE vote_usage() { return usage; };
	virtual bool vote_sampled() { return sampled; };
	
	virtual void render(const std::function<void(const std::vector<Vector3s>&, const std::vector<Vector3i>&, const Eigen::Quaternion<scalar>&, const Vector3s&, const scalar&)>&) const = 0;
	
	virtual void center(Vector3s& cent) const;
	
	inline bool is_root() const
	{
		return (parent == NULL);
	}
	
	DISTANCE_FIELD_TYPE type;
	DISTANCE_FIELD_USAGE usage;
	std::shared_ptr< DistanceField > parent;
	int group;
	int params_index;
	bool sampled;
};

struct DistanceFieldOperator : public DistanceField
{
	DistanceFieldOperator(DISTANCE_FIELD_TYPE type_, DISTANCE_FIELD_USAGE usage_, int group_, int params_index_, bool sampled_);
	
	virtual void advance(const scalar& dt);
	virtual scalar compute_phi_vel(const Vector3s& pos, Vector3s& vel) const;
	virtual scalar compute_phi(const Vector3s& pos) const;
	virtual bool local_bounding_box(Vector3s& low, Vector3s& high) const;
    virtual void resample_mesh(const scalar& dx, VectorXs& result, VectorXs& normals);
	virtual bool check_durations(const scalar& cur_time, const scalar& cur_vol, Vector3s& shooting_vel);
	virtual void apply_global_rotation(const Eigen::Quaternion<scalar>& rot);
	virtual void apply_local_rotation(const Eigen::Quaternion<scalar>& rot);
	virtual void apply_translation(const Vector3s& t);
	virtual int vote_param_indices();
	virtual DISTANCE_FIELD_USAGE vote_usage();
	virtual bool vote_sampled();

	virtual void render(const std::function<void(const std::vector<Vector3s>&, const std::vector<Vector3i>&, const Eigen::Quaternion<scalar>&, const Vector3s&, const scalar&)>&) const;
	
	std::vector< std::shared_ptr<DistanceField> > children;
};


struct DistanceFieldObject : public DistanceField
{
    DistanceFieldObject(const Vector3s& center_, const VectorXs& parameter_, DISTANCE_FIELD_TYPE type_, DISTANCE_FIELD_USAGE usage_, bool inside, const Vector3s& raxis, const scalar& rangle, int group_, int params_index_, bool sampled_, const std::vector< DF_SOURCE_DURATION >& durations_, const std::string& szfn = "", const std::string& szfn_cache = "");
	virtual void advance(const scalar& dt);
    virtual void resample_mesh(const scalar& dx, VectorXs& result, VectorXs& normals);
	virtual bool check_durations(const scalar& cur_time, const scalar& cur_vol, Vector3s& shooting_vel);
	virtual scalar compute_phi_vel(const Vector3s& pos, Vector3s& vel) const;
	virtual scalar compute_phi(const Vector3s& pos) const;
	virtual bool local_bounding_box(Vector3s& low, Vector3s& high) const;
	virtual void apply_global_rotation(const Eigen::Quaternion<scalar>& rot);
	virtual void apply_local_rotation(const Eigen::Quaternion<scalar>& rot);
	virtual void apply_translation(const Vector3s& t);
	
	virtual void render(const std::function<void(const std::vector<Vector3s>&, const std::vector<Vector3i>&, const Eigen::Quaternion<scalar>&, const Vector3s&, const scalar&)>&) const;
	
    void process_file_mesh(const std::string& szfn_cache);
    
	Vector3s center;
	VectorXs parameter;
	
	Eigen::Quaternion<scalar> rot;
	
	Vector3s future_center;
	Eigen::Quaternion<scalar> future_rot;
	
	Vector3s omega;
	Vector3s V;
	
	scalar sign;
	
	std::shared_ptr<SolidMesh> mesh;
    Array3s volume;
    Vector3s volume_origin;
	
	std::vector< DF_SOURCE_DURATION > durations;
};



#endif
