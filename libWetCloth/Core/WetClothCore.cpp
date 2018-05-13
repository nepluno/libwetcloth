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

#include "WetClothCore.h"
#include "TimingUtilities.h"
#include "MemUtilities.h"

WetClothCore::WetClothCore( const std::shared_ptr<TwoDScene>& scene, const std::shared_ptr<SceneStepper>& scene_stepper )
: m_scene(scene)
, m_scene_stepper(scene_stepper)
, m_current_step(0)
{
    timing_buffer.resize(15);
    timing_buffer.assign(15, 0.0);
    
    memset(&m_info, 0, sizeof(Info));
}

WetClothCore::~WetClothCore()
{}

const WetClothCore::Info& WetClothCore::getInfo() const
{
    return m_info;
}

const std::shared_ptr<SceneStepper>& WetClothCore::getSceneStepper() const
{
    return m_scene_stepper;
}

int WetClothCore::getCurrentTime() const
{
    return m_current_step;
}

const std::shared_ptr<TwoDScene>& WetClothCore::getScene() const
{
    return m_scene;
}

const std::vector<scalar>& WetClothCore::getTimingStatistics() const
{
    return timing_buffer;
}

void WetClothCore::stepSystem(const scalar &dt){
    assert( m_scene != NULL );
    assert( m_scene_stepper != NULL );
    
    VectorXs oldpos = m_scene->getX();
    VectorXs oldvel = m_scene->getV();
	
    const scalar max_elasto_vel = m_scene->getMaxVelocity();
    const scalar max_fluid_vel = m_scene->getMaxFluidVelocity();
    const scalar dx = m_scene->getCellSize();
    const scalar max_elasto_dt = std::min(dx / std::max(1e-63, max_elasto_vel) / 3.0, 1.0 / 30.0); // 1/6 CFLs
    const scalar max_fluid_dt = std::min(dx / std::max(1e-63, max_fluid_vel) * 3.0, 1.0 / 30.0); // 3 CFLs
    const scalar max_dt = std::min(max_elasto_dt, max_fluid_dt);
    
    const int num_substeps = std::max(1, (int) ceil(dt / max_dt));
    const scalar sub_dt = dt / (scalar) num_substeps;
	
	m_info.m_historical_max_vel = std::max(m_info.m_historical_max_vel, max_elasto_vel);
	m_info.m_historical_max_vel_fluid = std::max(m_info.m_historical_max_vel_fluid, max_fluid_vel);
    
    std::cout << "[step system max vel: (" << max_elasto_vel << " <" << m_info.m_historical_max_vel << ">, " << max_fluid_vel
    << " <" << m_info.m_historical_max_vel_fluid << ">) max dt: " << max_dt << " (" << max_elasto_dt << ", " << max_fluid_dt
    << "), # sub-step: (" << num_substeps << "), sub-dt: " << sub_dt << "]" << std::endl;
    
    for(int k = 0; k < num_substeps; ++k) {
		scalar cur_time = (scalar) m_current_step * dt + k * sub_dt;
        std::cout << "[(" << cur_time << " s) start substep: " << k << "/" << num_substeps << "]" << std::endl;

        scalar t0 = timingutils::seconds();
        scalar t1;
        
        m_scene->stepScript(sub_dt, cur_time);
        m_scene->applyScript(sub_dt);
		m_scene->sampleLiquidDistanceFields(cur_time + sub_dt);
        
        // Sample Grid around Particles
        m_scene->updateParticleBoundingBox();
        m_scene->rebucketizeParticles();
        m_scene->resampleNodes();
        t1 = timingutils::seconds();
        timing_buffer[0] += t1 - t0; // build Grid
        t0 = t1;
        
        // Apply Flow-Rule to Particles
        m_scene->postProcess(sub_dt);
        m_scene->updateManifoldOperators();
        m_scene->updateOrientation();
        m_scene->updateLiquidPhi(sub_dt);
        
        // Surface tension
        m_scene->updateIntersection();
        m_scene->updateStartState();
        m_scene->updateSolidPhi();
        m_scene->updateSolidWeights();
        
        // save velocity
        m_scene->saveParticleVelocity();
        t1 = timingutils::seconds();
        timing_buffer[1] += t1 - t0; // Plasticity, Solid Stress, and Compute Kernel Weights
        t0 = t1;
        
		m_scene->mapParticleNodesAPIC();
		m_scene->saveFluidVelocity();
		m_scene->mapParticleSaturationPsiNodes();
		m_scene->updatePorePressureNodes();
		t1 = timingutils::seconds();
		timing_buffer[2] += t1 - t0; // APIC Mapping
		t0 = t1;
		
		
		if(m_scene->useBiCGSTAB())
		{
			m_scene_stepper->stepVelocity( *m_scene, sub_dt );
            t1 = timingutils::seconds();
            timing_buffer[3] += t1 - t0; // Velocity Prediction
            t0 = t1;
			
			m_scene_stepper->solveBiCGSTAB( *m_scene, sub_dt );
			
			t1 = timingutils::seconds();
			timing_buffer[4] += t1 - t0;
			t0 = t1;
		} else {
			
			m_scene_stepper->stepVelocity( *m_scene, sub_dt );
            t1 = timingutils::seconds();
            timing_buffer[3] += t1 - t0; // Velocity Prediction
            t0 = t1;
            
            if(m_scene->getLiquidInfo().check_divergence) {
                m_info.m_initial_div_accu += m_scene_stepper->computeDivergence(*m_scene) / (scalar) num_substeps;
            }
            
			m_scene_stepper->projectFine( *m_scene, sub_dt );
			t1 = timingutils::seconds();
			timing_buffer[4] += t1 - t0; // Pressure Projection
			t0 = t1;
			
			if(m_scene->getLiquidInfo().solve_solid) {
				m_scene_stepper->applyPressureDragElasto(*m_scene, sub_dt);
                t1 = timingutils::seconds();
                timing_buffer[5] += t1 - t0; // Solve solid velocity
                t0 = t1;
            }
            
            if(m_scene->getLiquidInfo().check_divergence) {
                m_scene_stepper->pushFluidVelocity();
                m_scene_stepper->applyPressureDragFluid(*m_scene, sub_dt);
                scalar div = m_scene_stepper->computeDivergence(*m_scene) / (scalar) num_substeps;
                m_info.m_explicit_div_accu += div;
                m_scene_stepper->popFluidVelocity();
            }
                
            if(m_scene->getLiquidInfo().solve_solid) {
                m_scene_stepper->stepImplicitElasto( *m_scene, sub_dt );
                t1 = timingutils::seconds();
                timing_buffer[5] += t1 - t0; // Solve solid velocity
                t0 = t1;
            }
			
			m_scene_stepper->applyPressureDragFluid(*m_scene, sub_dt);
			t1 = timingutils::seconds();
			timing_buffer[6] += t1 - t0; // Solve Fluid velocity
			t0 = t1;
			
			m_scene_stepper->acceptVelocity(*m_scene);
            
            if(m_scene->getLiquidInfo().check_divergence) {
                m_info.m_implicit_div_accu += m_scene_stepper->computeDivergence(*m_scene) / (scalar) num_substeps;
            }
		}
		
		// Transfer Velocity Back to Particles
		m_scene->correctLiquidParticles(sub_dt);
		t1 = timingutils::seconds();
		timing_buffer[7] += t1 - t0; // Particle Correction
		t0 = t1;
		
		m_scene->mapNodeParticlesAPIC();
		t1 = timingutils::seconds();
		timing_buffer[8] += t1 - t0; // APIC Map Particle Back
		t0 = t1;
		
        // Advect Particles
        m_scene_stepper->advectScene( *m_scene, sub_dt );

        m_scene->solidProjection(dt);
        t1 = timingutils::seconds();
        timing_buffer[9] += t1 - t0; // Particle Advection
        t0 = t1;
		
        m_scene->distributeFluidElasto();
        t1 = timingutils::seconds();
        timing_buffer[10] += t1 - t0; // Liquid Capturing
        t0 = t1;
		
        m_scene->distributeElastoFluid();
        t1 = timingutils::seconds();
        timing_buffer[11] += t1 - t0; // Liquid Dripping
        t0 = t1;
        
        m_scene->updateVelocityDifference();
        
        m_scene->updateGaussAccel();
        
        m_scene_stepper->manifoldPropagate( *m_scene, sub_dt );
        t1 = timingutils::seconds();
        timing_buffer[12] += t1 - t0; // Solve Quasi-Static Equations
        t0 = t1;
        
        std::cout << "[terminate particles]" << std::endl;
        m_scene->terminateParticles();
        
        std::cout << "[update optimal volume]" << std::endl;
        m_scene->updateOptiVolume();
        
        std::cout << "[split particles]" << std::endl;
        m_scene->splitLiquidParticles();
        
        std::cout << "[merge particles]" << std::endl;
        m_scene->mergeLiquidParticles();
        t1 = timingutils::seconds();
        timing_buffer[13] += t1 - t0; // Merge & Split Particles
        t0 = t1;
        
        std::cout << "[update gauss system]" << std::endl;
        m_scene->updateGaussSystem(sub_dt);
        t1 = timingutils::seconds();
        timing_buffer[14] += t1 - t0; // update Deformation Gradient
        t0 = t1;
    }
    
	if(m_scene->getLiquidInfo().check_divergence) {
		scalar avg_explicit_div = m_info.m_explicit_div_accu / (scalar) (m_current_step + 1);
		scalar avg_implicit_div = m_info.m_implicit_div_accu / (scalar) (m_current_step + 1);
		scalar avg_initial_div = m_info.m_initial_div_accu / (scalar) (m_current_step + 1);
		std::cout << "Div Check, " << avg_initial_div << ", " << avg_explicit_div << ", " << avg_implicit_div << ", " << (fabs(avg_implicit_div - avg_explicit_div) / avg_initial_div) << std::endl;
	}
	
    size_t cur_usage = memutils::getCurrentRSS();
    
    m_info.m_mem_usage_accu += (scalar) cur_usage;
    m_info.m_num_particles_accu += (scalar) m_scene->getNumParticles();
    m_info.m_num_elements_accu += (scalar) m_scene->getNumGausses();
    m_info.m_num_fluid_particles_accu += (scalar) m_scene->getNumFluidParticles();

    // Check for obvious problems in the simulated scene
#ifdef DEBUG
    m_scene->checkConsistency();
#endif
	
	++m_current_step;
}



