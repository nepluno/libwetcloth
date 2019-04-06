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

/*
 * This is the main function where time stepping happens
 */
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
	
	// Start the possible sub-steps
    for(int k = 0; k < num_substeps; ++k) {
		scalar cur_time = (scalar) m_current_step * dt + k * sub_dt;
        std::cout << "[(" << cur_time << " s) start substep: " << k << "/" << num_substeps << "]" << std::endl;

        scalar t0 = timingutils::seconds();
        scalar t1;
		
		// Update Viscous Parameter for Elastic Rods
		m_scene->updateStrandParamViscosity(sub_dt);
		
		// Setup Scripting for Kinematic Objects
        m_scene->stepScript(sub_dt, cur_time);
        m_scene->applyScript(sub_dt);
		
		// Emit Liquid Particles for Liquid Sources
		m_scene->sampleLiquidDistanceFields(cur_time + sub_dt);
		
		// Remove Liquid Particles outside Simulation Domain (to save time)
		std::cout << "[terminate particles]" << std::endl;
		m_scene->terminateParticles();
		
		// Calculate the Optimal Volume of Liquid Particles
		std::cout << "[update optimal volume]" << std::endl;
		m_scene->updateOptiVolume();
		
		// Split the Liquid Particles if They are too Large
		std::cout << "[split particles]" << std::endl;
		m_scene->splitLiquidParticles();
		
		// Merge the Liquid Particles if They are too Small
		std::cout << "[merge particles]" << std::endl;
		m_scene->mergeLiquidParticles();
		t1 = timingutils::seconds();
		timing_buffer[0] += t1 - t0; // Merge & Split Particles
		t0 = t1;
		
        // Create Grid around Particles
        m_scene->updateParticleBoundingBox();
        m_scene->rebucketizeParticles();
        m_scene->resampleNodes();
        t1 = timingutils::seconds();
        timing_buffer[1] += t1 - t0; // build Grid
        t0 = t1;
        
        // Update Particle-Node Weight
        m_scene->computeWeights(sub_dt);
		
		// Update Solid Stress
		m_scene->computedEdFe();
		
        m_scene->updateManifoldOperators();

		// Update the Orientation Field
        m_scene->updateOrientation();
		
		// Update the Liquid Distance Field
        m_scene->updateLiquidPhi(sub_dt);
        
        // Compute Cohesion Force
        m_scene->updateIntersection();
        
        // Advect Surface Tension Force
        m_scene_stepper->advectSurfTension( *m_scene, dt );
		
		// Here's the precomputation of some forces lay
        m_scene->updateStartState();
		
		// Update the Distance Function for Kinematic Objects
        m_scene->updateSolidPhi();
		
		// Update the Weight on Grid (see [Batty et al. 2007] for details) for Kinematic Objects
        m_scene->updateSolidWeights();
        
        // Save Current Velocity
        m_scene->saveParticleVelocity();
		
        t1 = timingutils::seconds();
        timing_buffer[2] += t1 - t0; // Compute Weight, Solid Stress, and Distance Field (all above)
        t0 = t1;
		
		// Map the Liquid Particles and Elastic Vertices onto Grid
		m_scene->mapParticleNodesAPIC();
		
		// Save the Grid Velocity
		m_scene->saveFluidVelocity();
		
		// Update Saturation and Solid Volume Fraction on Grid
		m_scene->mapParticleSaturationPsiNodes();
		
		// Compute the Pore Pressure on Grid
		m_scene->updatePorePressureNodes();
		
		t1 = timingutils::seconds();
		timing_buffer[3] += t1 - t0; // APIC Mapping & Computing the Fields (all above)
		t0 = t1;
		
        // Explicitly Integrate the Elastic and Liquid Velocity
        m_scene_stepper->stepVelocity( *m_scene, sub_dt );
        t1 = timingutils::seconds();
        timing_buffer[4] += t1 - t0; // Velocity Prediction
        t0 = t1;
        
        // Check Divergence if Necessary
        if(m_scene->getLiquidInfo().check_divergence) {
            m_info.m_initial_div_accu += m_scene_stepper->computeDivergence(*m_scene) / (scalar) num_substeps;
        }
        
        // Do Pressure Projection for the Mixture
        m_scene_stepper->projectFine( *m_scene, sub_dt );
        t1 = timingutils::seconds();
        timing_buffer[5] += t1 - t0; // Pressure Projection
        t0 = t1;
        
        if(m_scene->getLiquidInfo().solve_solid) {
            // Apply Pressure Gradient to Solid
            m_scene_stepper->applyPressureDragElasto(*m_scene, sub_dt);
            t1 = timingutils::seconds();
            timing_buffer[6] += t1 - t0; // Timing the Pressure Gradient Application
            t0 = t1;
        }
        
        // Check Divergence if Necessary and Comparing with the Previously
        // Recorded Divergence to Measure the Error
        if(m_scene->getLiquidInfo().check_divergence) {
            m_scene_stepper->pushFluidVelocity();
            m_scene_stepper->applyPressureDragFluid(*m_scene, sub_dt);
            scalar div = m_scene_stepper->computeDivergence(*m_scene) / (scalar) num_substeps;
            m_info.m_explicit_div_accu += div;
            m_scene_stepper->popFluidVelocity();
        }
        
        if(m_scene->getLiquidInfo().solve_solid) {
            // Implicitly Integrate the Elastic Objects
            m_scene_stepper->stepImplicitElasto( *m_scene, sub_dt );
            t1 = timingutils::seconds();
            timing_buffer[6] += t1 - t0; // Solve solid velocity
            t0 = t1;
        }
        
        // Apply Pressure Gradient to Liquid
        m_scene_stepper->applyPressureDragFluid(*m_scene, sub_dt);
        t1 = timingutils::seconds();
        timing_buffer[7] += t1 - t0; // Solve Fluid velocity
        t0 = t1;
        
        // Update the Current Velocity with the Solved Ones
        m_scene_stepper->acceptVelocity(*m_scene);
        
        if(m_scene->getLiquidInfo().check_divergence) {
            m_info.m_implicit_div_accu += m_scene_stepper->computeDivergence(*m_scene) / (scalar) num_substeps;
        }
		
		// Kinematic Projection of the Liquid Velocity at the Boundary (as Fail-safe)
		m_scene->constrainLiquidVelocity();
		
		// Relax the Liquid Particles (see [Ando et al. 2011] for details)
		m_scene->correctLiquidParticles(sub_dt);
		t1 = timingutils::seconds();
		timing_buffer[8] += t1 - t0; // Particle Correction
		t0 = t1;
		
		// Transfer Velocity Back to Particles and Elastic Vertices
		m_scene->mapNodeParticlesAPIC();
		t1 = timingutils::seconds();
		timing_buffer[9] += t1 - t0; // APIC Map Particle Back
		t0 = t1;

		// Update the Multipliers applied on Geometric Stiffness
		// (refer to the supplemental material of [Fei et al. 2017] for details)
        m_scene->updateMultipliers( sub_dt );
		
        // Advection of Liquid Particles and Elastic Vertices
        m_scene_stepper->advectScene( *m_scene, sub_dt );

		// Kinematic Projection of the Elastic Vertices at the Boundary (as Fail-safe)
        m_scene->solidProjection( sub_dt );
        t1 = timingutils::seconds();
        timing_buffer[10] += t1 - t0; // Particle Advection
        t0 = t1;
		
		// Distribute the Liquid Volume onto Elastic Vertices (Capturing)
        m_scene->distributeFluidElasto(sub_dt);
        t1 = timingutils::seconds();
        timing_buffer[11] += t1 - t0; // Liquid Capturing
        t0 = t1;
		
		// Emit Liquid Particles for Overflowed Elastic Vertices (Dripping)
        m_scene->distributeElastoFluid();
        t1 = timingutils::seconds();
        timing_buffer[12] += t1 - t0; // Liquid Dripping
        t0 = t1;
		
		// Update the Velocity Displacement
        m_scene->updateVelocityDifference();
		
		// Update the Acceleration of Liquid on Elastic Vertices
        m_scene->updateGaussAccel();
		
		// Solve the Quasi-Static Equation on Elastic Vertices
        m_scene_stepper->manifoldPropagate( *m_scene, sub_dt );
        t1 = timingutils::seconds();
        timing_buffer[13] += t1 - t0; // Solve Quasi-Static Equations
        t0 = t1;
		
		// Update the Variables on Elements
		// We denote elements as 'Gauss' since they are computed at the Gaussian Quadrature Point (1-Point).
        std::cout << "[update gauss system and plasticity]" << std::endl;
        m_scene->updateGaussSystem(sub_dt);
		m_scene->updatePlasticity(sub_dt);
        t1 = timingutils::seconds();
        timing_buffer[14] += t1 - t0; // update Deformation Gradient
        t0 = t1;
    }
	
	// Summarize Divergence if Necessary
	if(m_scene->getLiquidInfo().check_divergence) {
		scalar avg_explicit_div = m_info.m_explicit_div_accu / (scalar) (m_current_step + 1);
		scalar avg_implicit_div = m_info.m_implicit_div_accu / (scalar) (m_current_step + 1);
		scalar avg_initial_div = m_info.m_initial_div_accu / (scalar) (m_current_step + 1);
		std::cout << "Div Check, " << avg_initial_div << ", " << avg_explicit_div << ", " << avg_implicit_div << ", " << (fabs(avg_implicit_div - avg_explicit_div) / avg_initial_div) << std::endl;
	}
	
	// Summarize Memory Usage
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



