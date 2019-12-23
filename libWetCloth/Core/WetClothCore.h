//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __WET_CLOTH_CORE_H__
#define __WET_CLOTH_CORE_H__

#include "TwoDScene.h"
#include "SceneStepper.h"

class WetClothCore
{
public:
    struct Info {
        scalar m_mem_usage_accu;
        scalar m_num_particles_accu;
        scalar m_num_fluid_particles_accu;
        scalar m_num_elements_accu;
        scalar m_initial_div_accu;
        scalar m_explicit_div_accu;
        scalar m_implicit_div_accu;
        scalar m_historical_max_vel;
        scalar m_historical_max_vel_fluid;
    };
    
    WetClothCore( const std::shared_ptr<TwoDScene>& scene, const std::shared_ptr<SceneStepper>& scene_stepper );

    virtual ~WetClothCore();
    /////////////////////////////////////////////////////////////////////////////
    // Simulation Control Functions
    
    virtual void stepSystem( const scalar& dt );
    
    virtual const std::vector<scalar>& getTimingStatistics() const;
    virtual const std::shared_ptr<TwoDScene>& getScene() const;
    virtual const std::shared_ptr<SceneStepper>& getSceneStepper() const;
    virtual const Info& getInfo() const;
    
    virtual int getCurrentTime() const;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	std::shared_ptr<TwoDScene> m_scene;
    std::shared_ptr<SceneStepper> m_scene_stepper;
    
	int m_current_step;
    
    std::vector<scalar> timing_buffer;
    
    Info m_info;
};

#endif
