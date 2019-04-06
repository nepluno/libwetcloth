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
