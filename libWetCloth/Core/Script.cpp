//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "Script.h"
#include "MathUtilities.h"
#include "TwoDScene.h"
#include "DistanceFields.h"

using namespace mathutils;

scalar Script::getNextVelocity( const scalar& dt, const scalar& current_time )
{
	if(func == Script::COSINE) {
		return cosine_ease_function(current_time + dt, start, end, start + ease_start, end - ease_end, amplitude, frequency);
	} else if(func == Script::CUBIC) {
		return cubic_ease_function(current_time + dt, start, end, start + ease_start, end - ease_end, v(3));
	} else if(func == Script::WENO) {
		const scalar vel = weno_ease_function(current_time + dt, dt, start, end, base_dt, base_pos, base_vertices);
		base_pos += vel * dt;
		return vel;
	}
	
	return 0.0;
}

void Script::stepScript( const scalar& dt, const scalar& current_time )
{
	if(!m_scene) return;
	
	if(current_time >= start && current_time + dt <= end)
	{
		switch (type) {
			case Script::TRANSLATE:
			{
				Vector3s vec(v(0), v(1), v(2));
				scalar dx = vec.norm();
				v(3) = dx;
				
				Vector3s& trans = m_scene->getGroupTranslation( group_index );
				m_scene->getPrevGroupTranslation( group_index ) = trans;
				if(dx > 0.0) {
					VectorXs nv = vec / dx;
					
					scalar vel = getNextVelocity( dt, current_time );
					trans += nv * vel * dt;
					m_scene->getGroupDistanceField( group_index )->apply_translation(nv * vel * dt);
				}

				break;
			}
			case Script::ROTATE:
			{
				scalar vel = getNextVelocity( dt, current_time );
				scalar da = vel * dt;
				Vector3s axis(v(0), v(1), v(2));
				
				Eigen::AngleAxis<scalar> rot(da, axis);
				Eigen::Quaternion<scalar> qrot(rot);
				
				Eigen::Quaternion<scalar>& prot = m_scene->getGroupRotation( group_index );
				m_scene->getPrevGroupRotation( group_index ) = prot;
                
                if(transform_with_origin) {
                    Vector3s& trans = m_scene->getGroupTranslation( group_index );
                    m_scene->getPrevGroupTranslation( group_index ) = trans;
                }
                
				if(transform_global) {
					prot = qrot * prot;
					m_scene->getGroupDistanceField( group_index )->apply_global_rotation(qrot);
				} else {
					prot *= qrot;
					m_scene->getGroupDistanceField( group_index )->apply_local_rotation(qrot);
				}
                
                if(transform_with_origin) {
                    const Eigen::Quaternion<scalar> q_diff = prot * m_scene->getPrevGroupRotation( group_index ).inverse();
                    Vector3s& trans = m_scene->getGroupTranslation( group_index );
                    Vector3s dir_vec = trans - origin;
                    
                    Vector3s rot_dir_vec = q_diff * dir_vec;
                    Vector3s vec = rot_dir_vec - dir_vec;
                    trans += vec;
                    m_scene->getGroupDistanceField( group_index )->apply_translation(vec);
                }
				break;
			}
			default:
				std::cout << "UNIMPLEMENTED SCRIPT TYPE [" << type << "]!" << std::endl;
				break;
		}
	}
}
