//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __SCRIPT_H__
#define __SCRIPT_H__

#include <Eigen/Core>
#include <memory>
#include "MathDefs.h"

class TwoDScene;

struct Script
{
	enum TYPE
	{
		ROTATE,
		TRANSLATE,
		
		TYPE_COUNT
	};
	
	enum FUNC
	{
		CUBIC,
		COSINE,
		WENO,
		
		FUNC_COUNT
	};
	
	TYPE type;
	FUNC func;
	int group_index;
	Vector4s v;
	Vector3s origin;
	scalar start;
	scalar end;
	scalar ease_start;
	scalar ease_end;
	scalar amplitude;
	scalar frequency;
	scalar base_dt;
	scalar base_pos;
	
	std::shared_ptr<TwoDScene> m_scene;

	bool transform_global;
    bool transform_with_origin;
	
	std::vector<scalar> base_vertices;
	
	scalar getNextVelocity( const scalar& dt, const scalar& current_time );
	
	void stepScript( const scalar& dt, const scalar& current_time );
	
};


#endif
