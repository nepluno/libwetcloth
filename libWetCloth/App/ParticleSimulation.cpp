//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ParticleSimulation.h"
#include "TimingUtilities.h"
#include "MemUtilities.h"

#ifdef RENDER_ENABLED
#include <AntTweakBar.h>
#endif

const static char* g_timing_labels[] = {
	"Sample, Merge and Split Particles",
    "Build Sparse MAC Grid",
    "Compute Weight, Solid Stress, and Distance Field",
    "APIC Particle/Vertex-to-Grid Transfer",
    "Force Integration and Velocity Prediction",
    "Solve Poisson Equation",
    "Solve Solid Velocity",
    "Solve Liquid Velocity",
    "Particle Correction",
    "APIC Grid-to-Particle/Vertex Transfer",
    "Particle/Vertex Advection",
    "Liquid Capturing",
    "Liquid Dripping",
    "Solve Quasi-Static Equation",
    "Update Deformation Gradient and Plasticity"
};
const static int num_timing_labels = sizeof(g_timing_labels) / sizeof(const char*);

ParticleSimulation::ParticleSimulation( const std::shared_ptr<TwoDScene>& scene, const std::shared_ptr<SceneStepper>& scene_stepper, const std::shared_ptr<TwoDSceneRenderer>& scene_renderer )
: m_core(std::make_shared<WetClothCore>( scene, scene_stepper ))
, m_scene_renderer(scene_renderer)
, m_display_controller( std::make_shared<TwoDimensionalDisplayController>( 1920, 1080 ) )
{
}

ParticleSimulation::~ParticleSimulation(){
}

int ParticleSimulation::currentCameraIndex() const
{
    return m_display_controller->currentCameraIndex();
}

void ParticleSimulation::stepSystem(const scalar &dt)
{
    m_core->stepSystem(dt);
    
    const std::vector<scalar>& timing_buffer = m_core->getTimingStatistics();
    
    scalar total_time = 0.0;
    for(scalar t : timing_buffer) total_time += t;
    
    std::cout << "---------------------------------" << std::endl;
    for(int i = 0; i < num_timing_labels; ++i) {
        scalar avg_time = (timing_buffer[i] / (scalar) (m_core->getCurrentTime() + 1));
        scalar prop = timing_buffer[i] / total_time * 100.0;
        std::cout << g_timing_labels[i] << ", " << avg_time << ", " << prop << "%" << std::endl;
    }
    
    const scalar divisor = (scalar) (m_core->getCurrentTime() + 1);
    
    std::cout << "---------------------------------" << std::endl;
    std::cout << "Total Time (per Frame), " << total_time << ", " << (total_time / divisor) << std::endl;
    scalar part_fluid_vol = m_core->getScene()->totalFluidVolumeParticles();
    scalar vert_fluid_vol = m_core->getScene()->totalFluidVolumeSoftElasto();
    std::cout << "Liquid Vol, " << part_fluid_vol << ", " << vert_fluid_vol << ", " << (part_fluid_vol + vert_fluid_vol) << std::endl;
    
    
    int peak_idx = 0;
    int cur_idx = 0;
    
    const char* mem_units[] = {
        "B", "KB", "MB", "GB", "TB", "PB"
    };
    
    size_t peak_usage = memutils::getPeakRSS();
    scalar peak_mem = (scalar) peak_usage;
    while(peak_mem > 1024.0 && peak_idx < (int)(sizeof(mem_units) / sizeof(char*))) {peak_mem /= 1024.0; peak_idx++;}
    
    const WetClothCore::Info& info = m_core->getInfo();
    
    scalar avg_mem = info.m_mem_usage_accu / divisor;
    while(avg_mem > 1024.0 && cur_idx < (int)(sizeof(mem_units) / sizeof(char*))) {avg_mem /= 1024.0; cur_idx++;}
    
    std::cout << "Particles (Avg.), " << m_core->getScene()->getNumParticles() << ", " << (info.m_num_particles_accu / divisor) << ", Fluid (Avg.), " << m_core->getScene()->getNumFluidParticles() << ", " << (info.m_num_fluid_particles_accu / divisor) << ", Elements, " << m_core->getScene()->getNumGausses() << ", " << (info.m_num_elements_accu / divisor) << std::endl;
    std::cout << "Peak Mem Usage, " << peak_mem << mem_units[peak_idx] << ", Avg Mem Usage, " << avg_mem << mem_units[cur_idx] << std::endl;
    
    std::cout << "---------------------------------" << std::endl;
}

void ParticleSimulation::initializeOpenGLRenderer(){
#ifdef RENDER_ENABLED
    TwBar* bar = TwNewBar("Control Panel");
    TwDefine(" GLOBAL help='A Multi-Scale Model for Simulating Liquid-Fabric Interactions' ");
    TwDefine(" 'Control Panel' size='400 600' color='83 37 0' position='5 70'");
    
    RenderInfo& render = m_scene_renderer->getRenderInfo();
    TwAddVarRW(bar, "particles", TW_TYPE_BOOLCPP, &render.render_particles, " group='visualization'");
    TwAddVarRW(bar, "vertices", TW_TYPE_BOOLCPP, &render.render_vertices, " group='visualization'");
    TwAddVarRW(bar, "elements", TW_TYPE_BOOLCPP, &render.render_gauss, " group='visualization'");
    
    TwAddSeparator(bar, NULL, " group='visualization'");
    TwAddVarRW(bar, "particle velocity", TW_TYPE_BOOLCPP, &render.render_particle_velocity, " group='visualization'");
    TwAddVarRW(bar, "vertex velocity", TW_TYPE_BOOLCPP, &render.render_vertice_velocity, " group='visualization'");
    TwAddVarRW(bar, "element velocity", TW_TYPE_BOOLCPP, &render.render_gauss_velocity, " group='visualization'");
    TwAddVarRW(bar, "velocity length", TW_TYPE_DOUBLE, &render.render_velocity_length, " min=0.0 step=0.1 group='visualization'");
    TwAddVarRW(bar, "deformation gradient length", TW_TYPE_DOUBLE, &render.render_deformation_gradient_length, " min=0.0 step=0.01 group='visualization'");
    TwAddVarRW(bar, "cohesion", TW_TYPE_BOOLCPP, &render.render_cohesion, " group='visualization'");
    
    TwAddSeparator(bar, NULL, " group='visualization'");
    TwAddVarRW(bar, "cloth", TW_TYPE_BOOLCPP, &render.render_cloth, " group='visualization'");
    TwAddVarRW(bar, "yarn", TW_TYPE_BOOLCPP, &render.render_yarn, " group='visualization'");
    TwAddVarRW(bar, "levelset", TW_TYPE_BOOLCPP, &render.render_levelset, " group='visualization'");
    TwAddVarRW(bar, "spring", TW_TYPE_BOOLCPP, &render.render_spring, " group='visualization'");
    
    TwAddSeparator(bar, NULL, " group='visualization'");
    TwAddVarRW(bar, "buckets", TW_TYPE_BOOLCPP, &render.render_buckets, " group='visualization'");
    {
        TwEnumVal EV[] = { {RenderInfo::RI_NV_NONE, "hide"}, {RenderInfo::RI_NV_CONSTANT, "show"}, {RenderInfo::RI_NV_SOLID_PHI, "solid SDF"} };
        TwType Type = TwDefineEnum("NVType", EV, 3);
        TwAddVarRW(bar, "nodes", Type, &render.render_nodes, " group='visualization'");
    }
    
    {
        TwEnumVal EV[] = { {RenderInfo::RI_FCV_NONE, "hide"}, {RenderInfo::RI_FCV_CONSTANT, "show"}, {RenderInfo::RI_FCV_SOLID_VOL, "volume fraction"}, {RenderInfo::RI_FCV_LIQUID_VOL, "saturation"} };
        TwType Type = TwDefineEnum("FCVType", EV, 4);
        TwAddVarRW(bar, "face centers", Type, &render.render_face_centers, " group='visualization'");
    }
    
    {
        TwEnumVal EV[] = { {RenderInfo::RI_ECV_NONE, "hide"}, {RenderInfo::RI_ECV_CONSTANT, "show"} };
        TwType Type = TwDefineEnum("ECVType", EV, 2);
        TwAddVarRW(bar, "edge centers", Type, &render.render_edge_centers, " group='visualization'");
    }
    
    {
        TwEnumVal EV[] = { {RenderInfo::RI_CCV_NONE, "hide"}, {RenderInfo::RI_CCV_CONSTANT, "show"}, {RenderInfo::RI_CCV_LIQUID_PHI, "liquid SDF"} };
        TwType Type = TwDefineEnum("CCVType", EV, 3);
        TwAddVarRW(bar, "cell centers", Type, &render.render_cell_centers, " group='visualization'");
    }
    
    LiquidInfo& info = m_core->getScene()->getLiquidInfo();
    
    TwAddVarRO(bar, "pore radius", TW_TYPE_DOUBLE, &info.pore_radius, " help='Pore radius (cm)' group='static parameters'");
    TwAddVarRO(bar, "fiber diameter", TW_TYPE_DOUBLE, &info.yarn_diameter, " help='Yarn (mesh-based) / fiber (yarn-based) diameter (cm)' group='static parameters'");
    TwAddVarRO(bar, "rest volume fraction", TW_TYPE_DOUBLE, &info.rest_volume_fraction, " help='Volume fraction at rest' group='static parameters'");
    
    TwAddVarRW(bar, "liquid density", TW_TYPE_DOUBLE, &info.liquid_density, " min=0.1 max=10.0 step=0.1 help='Liquid density  (g/cm^3)' group='dynamic parameters'");
    TwAddVarRW(bar, "air density", TW_TYPE_DOUBLE, &info.air_density, " min=0.0 max=0.1 step=1e-3 help='Air density  (g/cm^3)' group='dynamic parameters'");
    TwAddVarRW(bar, "surface tension", TW_TYPE_DOUBLE, &info.surf_tension_coeff, " min=0.0 max=400.0 step=0.1 help='Surface tension coefficient (dyn/cm)' group='dynamic parameters'");
    TwAddVarRW(bar, "liquid viscosity", TW_TYPE_DOUBLE, &info.viscosity, " min=0.0 max=100.0 step=0.0001 help='Liquid viscosity (dyn·s/cm^2)' group='dynamic parameters'");
    TwAddVarRW(bar, "air viscosity", TW_TYPE_DOUBLE, &info.air_viscosity, " min=0.0 max=100.0 step=0.0001 help='Air viscosity (dyn·s/cm^2)' group='dynamic parameters'");
    TwAddVarRW(bar, "contact angle", TW_TYPE_DOUBLE, &info.rest_contact_angle, " min=0.0 max=1.8325957146 step=0.0001 help='Liquid-solid contact angle (rad)' group='dynamic parameters'");
    TwAddVarRW(bar, "drag nonlinearity", TW_TYPE_DOUBLE, &info.yazdchi_power, " min=0.0 max=3.0 step=0.1 help='Drag nonlinearity' group='dynamic parameters'");
    TwAddVarRW(bar, "power of fraction", TW_TYPE_DOUBLE, &info.lambda, " min=0.0 max=3.0 step=0.1 help='Volume Fraction Power' group='dynamic parameters'");
    TwAddVarRW(bar, "cohesion multiplier", TW_TYPE_DOUBLE, &info.cohesion_coeff, " min=0.0 max=1.0 step=0.001 help='Strength of cohesion' group='dynamic parameters'");
    TwAddVarRW(bar, "correction range", TW_TYPE_DOUBLE, &info.correction_multiplier, " min=0.0 max=10.0 step=0.1 help='Range of liquid particle correction' group='dynamic parameters'");
    TwAddVarRW(bar, "correction strength", TW_TYPE_DOUBLE, &info.correction_strength, " min=0.0 max=10.0 step=0.1 help='Strength of liquid particle correction' group='dynamic parameters'");
    TwAddVarRW(bar, "steps per correction", TW_TYPE_INT32, &info.correction_step, " min=1 max=1000 step=1 help='# steps between each time of correction' group='dynamic parameters'");
    TwAddVarRW(bar, "passes of curvature smoothing", TW_TYPE_INT32, &info.surf_tension_smoothing_step, " min=0 max=1000 step=1 help='# passes of Laplacian smoothing of the curvature before computing surface tension' group='dynamic parameters'");
    TwAddVarRW(bar, "Liquid 1st-order lossless", TW_TYPE_DOUBLE, &info.flip_coeff, " min=0.0 max=1.0 step=0.0001 help='Lossless of liquid 1st-order movement' group='dynamic parameters'");
    TwAddVarRW(bar, "Solid rotation/shearing lossless", TW_TYPE_DOUBLE, &info.elasto_flip_asym_coeff, " min=0.0 max=1.0 step=0.0001 help='Lossless of solid rotation/shearing' group='dynamic parameters'");
    TwAddVarRW(bar, "Solid stretching lossless", TW_TYPE_DOUBLE, &info.elasto_flip_coeff, " min=0.0 max=1.0 step=0.0001 help='Lossless of solid stretching' group='dynamic parameters'");
    TwAddVarRW(bar, "Solid advection lossless", TW_TYPE_DOUBLE, &info.elasto_advect_coeff, " min=0.0 max=1.0 step=0.0001 help='Lossless of solid advection' group='dynamic parameters'");
    TwAddVarRW(bar, "Levelset modulus", TW_TYPE_DOUBLE, &info.levelset_young_modulus, " min=0.0 max=1e8 step=1e4 help='Levelset Youngs modulus (barye)' group='dynamic parameters'");
    {
        TwEnumVal bendingEV[] = { {0, "tan(angle)"}, {1, "sin(angle)"}, {2, "angle"} };
        TwType bendingType = TwDefineEnum("bendingType", bendingEV, 3);
        TwAddVarRW(bar, "Bending energy", bendingType, &info.bending_scheme, " help='Change the way to calculate bending energy of cloth' group='dynamic parameters'");
    }
    
    TwAddVarRW(bar, "Use surface tension", TW_TYPE_BOOLCPP, &info.use_surf_tension, " help='Use/Not use surface tension force' group='features'");
    TwAddVarRW(bar, "Use cohesion", TW_TYPE_BOOLCPP, &info.use_cohesion, " help='Use/Not use cohesion force' group='features'");
    TwAddVarRW(bar, "Levelset cohesion", TW_TYPE_BOOLCPP, &info.solid_cohesion, " help='Calculate cohesion force between cloth/yarn and levelset' group='features'");
    TwAddVarRW(bar, "Cloth/Yarn cohesion", TW_TYPE_BOOLCPP, &info.soft_cohesion, " help='Calculate self-cohesion of cloth/yarn' group='features'");
    TwAddVarRW(bar, "Solve cloth/yarn", TW_TYPE_BOOLCPP, &info.solve_solid, " help='Calculate the dynamics of cloth/yarn' group='features'");
    TwAddVarRW(bar, "Use nonlinear drag", TW_TYPE_BOOLCPP, &info.use_nonlinear_drag, " help='Use nonlinear drag force' group='features'");
    TwAddVarRW(bar, "Use levelset collision", TW_TYPE_BOOLCPP, &info.use_levelset_force, " help='Apply penalty force between cloth/yarn and levelset' group='features'");
    TwAddVarRW(bar, "Pressure gradient on manifold", TW_TYPE_BOOLCPP, &info.apply_pressure_manifold, " help='Apply pressure gradient from liquid when solving flows on manifold' group='features'");
    TwAddVarRW(bar, "Use twisting force", TW_TYPE_BOOLCPP, &info.use_twist, " help='Apply twisting force' group='features'");
    TwAddVarRW(bar, "Use BiCGSTAB", TW_TYPE_BOOLCPP, &info.use_bicgstab, " help='Solve dynamics with BiCGSTAB in a brute-force style' group='features'");
    TwAddVarRW(bar, "Use AMGPCG for cloth/yarn", TW_TYPE_BOOLCPP, &info.use_amgpcg_solid, " help='Solve cloth/yarn dynamics with AMGPCG solver' group='features'");
    TwAddVarRW(bar, "Use PCR for cloth/yarn", TW_TYPE_BOOLCPP, &info.use_pcr, " help='Solve cloth/yarn dynamics with preconditioned conjugate residual solver (turn off to use conjugate gradient)' group='features'");
    TwAddVarRW(bar, "Pore pressure deforms cloth/yarn", TW_TYPE_BOOLCPP, &info.apply_pore_pressure_solid, " help='Cloth/yarn dynamics are affected by pore pressure' group='features'");
    TwAddVarRW(bar, "Propagate cloth/yarn velocity", TW_TYPE_BOOLCPP, &info.propagate_solid_velocity, " help='Propagate cloth/yarn velocity when solving flows on manifold' group='features'");
    TwAddVarRW(bar, "Check divergence", TW_TYPE_BOOLCPP, &info.check_divergence, " help='Check the divergence after pressure projection' group='features'");
    TwAddVarRW(bar, "Update volume fraction", TW_TYPE_BOOLCPP, &info.use_varying_fraction, " help='Update the volume fraction of cloth/yarn after solving for their dynamics' group='features'");
    TwAddVarRW(bar, "Compute viscosity", TW_TYPE_BOOLCPP, &info.use_varying_fraction, " help='Solve the viscosity term for liquid' group='features'");
    TwAddVarRW(bar, "Implicit viscosity", TW_TYPE_BOOLCPP, &info.implicit_viscosity, " help='Use implicit viscosity computation' group='features'");
    TwAddVarRW(bar, "Implicit cloth/yarn for liquid drag", TW_TYPE_BOOLCPP, &info.drag_by_future_solid, " help='Use implicitly (turn off to use explicitly integrated cloth/yarn velocity) integrated cloth/yarn velocity to compute liquid drag force' group='features'");
    TwAddVarRW(bar, "Air drag", TW_TYPE_BOOLCPP, &info.drag_by_future_solid, " help='Turn on to apply air drag force to the cloth/yarn' group='features'");
    TwAddVarRO(bar, "Init with nonuniform volume fraction", TW_TYPE_BOOLCPP, &info.init_nonuniform_fraction, " help='Must be turned on when initializing the cloth/yarn with nonuniform volume fraction' group='features'");
#endif
}


void ParticleSimulation::renderSceneOpenGL(const scalar& dt){
    assert( m_scene_renderer != NULL );
    m_scene_renderer->renderParticleSimulation(*m_core->getScene(), dt);
    
}

void ParticleSimulation::updateOpenGLRendererState(){
    assert( m_scene_renderer != NULL );
    m_scene_renderer->updateParticleSimulationState(*m_core->getScene());
}

void ParticleSimulation::computeCameraCenter(renderingutils::Viewport &view){
    const VectorXs& x = m_core->getScene()->getX();
    
    // Compute the bounds on all particle positions
    scalar max_x = -std::numeric_limits<scalar>::infinity();
    scalar min_x =  std::numeric_limits<scalar>::infinity();
    scalar max_y = -std::numeric_limits<scalar>::infinity();
    scalar min_y =  std::numeric_limits<scalar>::infinity();
	scalar max_z = -std::numeric_limits<scalar>::infinity();
	scalar min_z =  std::numeric_limits<scalar>::infinity();
    for( int i = 0; i < m_core->getScene()->getNumParticles(); ++i )
    {
        if( x(4*i) > max_x )   max_x = x(4*i);
        if( x(4*i) < min_x )   min_x = x(4*i);
        if( x(4*i+1) > max_y ) max_y = x(4*i+1);
        if( x(4*i+1) < min_y ) min_y = x(4*i+1);
		if( x(4*i+2) > max_z ) max_z = x(4*i+2);
		if( x(4*i+2) < min_z ) min_z = x(4*i+2);
    }
	
	const auto& gdf = m_core->getScene()->getGroupDistanceField();
	for(const auto& df : gdf)
	{
		Vector3s low, high;
		if(df->local_bounding_box(low, high)) continue;
		max_x = std::max(high(0), max_x);
		max_y = std::max(high(1), max_y);
		max_z = std::max(high(2), max_z);
		
		min_x = std::min(low(0), min_x);
		min_y = std::min(low(1), min_y);
		min_z = std::min(low(2), min_z);
	}
	
    // Set center of view to center of bounding box
    view.cx = 0.5*(max_x+min_x);
    view.cy = 0.5*(max_y+min_y);
	view.cz = 0.5*(max_z+min_z);
    
    // Set the zoom such that all particles are in view
    view.rx = 0.5*(max_x-min_x);
    if( view.rx == 0.0 ) view.rx = 1.0;
    view.ry = 0.5*(max_y-min_y);
    if( view.ry == 0.0 ) view.ry = 1.0;
	view.rz = 0.5*(max_z-min_z);
	if( view.rz == 0.0 ) view.rz = 1.0;
}

void ParticleSimulation::readPos( const std::string& fn_pos )
{
    std::ifstream ifs(fn_pos, std::ios::binary);
    m_scene_serializer.loadPosOnly(*m_core->getScene(), ifs);
    ifs.close();
}

void ParticleSimulation::serializePositionOnly( const std::string& fn_pos )
{
    m_scene_serializer.serializePositionOnly(*m_core->getScene(), fn_pos);
}

void ParticleSimulation::serializeScene(const std::string& fn_clothes,
                                        const std::string& fn_hairs,
                                        const std::string& fn_fluid,
                                        const std::string& fn_internal_boundaries,
                                        const std::string& fn_external_boundaries,
                                        const std::string& fn_spring){
    m_scene_serializer.serializeScene( *m_core->getScene(), fn_clothes, fn_hairs, fn_fluid, fn_internal_boundaries, fn_external_boundaries, fn_spring );
}

void ParticleSimulation::centerCamera(bool b_reshape)
{
	renderingutils::Viewport view;
	
	computeCameraCenter(view);
	scalar ratio;
	
	ratio = ((scalar)m_display_controller->getWindowHeight())/((scalar)m_display_controller->getWindowWidth());
	
	view.size = 1.2*std::max(ratio*view.rx,view.ry);
	
	m_display_controller->setCenterX(view.cx);
	m_display_controller->setCenterY(view.cy);
	m_display_controller->setCenterZ(view.cz);
	m_display_controller->setScaleFactor(view.size);
	m_display_controller->initCamera(view);
	
	if(b_reshape) m_display_controller->reshape(m_display_controller->getWindowWidth(),m_display_controller->getWindowHeight());
}

std::string ParticleSimulation::getSolverName(){
    return m_core->getSceneStepper()->getName();
}

void ParticleSimulation::keyboard( unsigned char key, int x, int y )
{
	m_display_controller->keyboard(key,x,y);
}


void ParticleSimulation::reshape( int w, int h )
{
	m_display_controller->reshape(w,h);
}


void ParticleSimulation::special( int key, int x, int y )
{
	m_display_controller->special(key, x, y);
}


void ParticleSimulation::mouse( int button, int state, int x, int y )
{
	m_display_controller->mouse(button,state,x,y);
}

void ParticleSimulation::motion( int x, int y )
{
	m_display_controller->motion(x,y);
}


int ParticleSimulation::getWindowWidth() const
{
	return m_display_controller->getWindowWidth();
}


int ParticleSimulation::getWindowHeight() const
{
	return m_display_controller->getWindowHeight();
}


void ParticleSimulation::setWindowWidth(int w)
{
	m_display_controller->setWindowWidth(w);
}

void ParticleSimulation::setWindowHeight(int h)
{
	m_display_controller->setWindowHeight(h);
}

const std::shared_ptr<TwoDimensionalDisplayController>& ParticleSimulation::getDC() const
{
	return m_display_controller;
}

void ParticleSimulation::setCenterX( double x )
{
	m_display_controller->setCenterX(x);
}

void ParticleSimulation::setCenterY( double y )
{
	m_display_controller->setCenterY(y);
}

void ParticleSimulation::setCenterZ( double z )
{
	m_display_controller->setCenterZ(z);
}

void ParticleSimulation::setScaleFactor( double scale )
{
	m_display_controller->setScaleFactor(scale);
}

void ParticleSimulation::setCamera( const Camera& cam )
{
	m_display_controller->getCamera() = cam;
	
	m_display_controller->setCenterX(cam.center_(0));
	m_display_controller->setCenterY(cam.center_(1));
	m_display_controller->setCenterZ(cam.center_(2));
	
	scalar size = 1.2 * cam.radius_;
	m_display_controller->setScaleFactor(size);
}

void ParticleSimulation::setView( const renderingutils::Viewport& view )
{
	m_display_controller->setCenterX(view.cx);
	m_display_controller->setCenterY(view.cy);
	m_display_controller->setCenterZ(view.cz);
	m_display_controller->setScaleFactor(view.size);
	m_display_controller->initCamera(view);

}

void ParticleSimulation::printDDA()
{
	m_core->getScene()->computeDDA();
}

void ParticleSimulation::finalInit()
{
    m_scene_serializer.initializeFaceLoops(*m_core->getScene());
}

const LiquidInfo& ParticleSimulation::getLiquidInfo()
{
    return m_core->getScene()->getLiquidInfo();
}

