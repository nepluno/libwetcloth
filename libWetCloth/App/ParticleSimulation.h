//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __PARTICLE_SIMULATION_H__
#define __PARTICLE_SIMULATION_H__

#include "TwoDScene.h"
#include "TwoDSceneRenderer.h"
#include "SceneStepper.h"
#include "TwoDSceneSerializer.h"
#include "TwoDimensionalDisplayController.h"
#include "WetClothCore.h"

// TODO: Move code out of header!
extern bool g_rendering_enabled;

class ParticleSimulation
{
public:
	
    ParticleSimulation( const std::shared_ptr<TwoDScene>& scene, const std::shared_ptr<SceneStepper>& scene_stepper, const std::shared_ptr<TwoDSceneRenderer>& scene_renderer);
    
    virtual ~ParticleSimulation();
    /////////////////////////////////////////////////////////////////////////////
    // Simulation Control Functions
    
    void stepSystem( const scalar& dt );
    
    /////////////////////////////////////////////////////////////////////////////
    // Rendering Functions
    
    void initializeOpenGLRenderer();
    void renderSceneOpenGL(const scalar& dt);
    void updateOpenGLRendererState();
    void computeCameraCenter(renderingutils::Viewport& view);


    /////////////////////////////////////////////////////////////////////////////
    // Serialization Functions
    void serializeScene( const std::string& fn_clothes,
                        const std::string& fn_hairs,
                        const std::string& fn_fluid,
                        const std::string& fn_internal_boundaries,
                        const std::string& fn_external_boundaries,
                        const std::string& fn_spring);
    
    void serializePositionOnly( const std::string& fn_pos );
    
    void readPos( const std::string& fn_pos );
    /////////////////////////////////////////////////////////////////////////////
    // Status Functions
    
	std::string getSolverName();
	
	void centerCamera(bool b_reshape = true);
	void keyboard( unsigned char key, int x, int y );
	void reshape( int w, int h );
	
	void special( int key, int x, int y );
	
	void mouse( int button, int state, int x, int y );
	
	void translateView( double dx, double dy );
	
	void zoomView( double dx, double dy );
	
	void motion( int x, int y );
	
	int getWindowWidth() const;
	int getWindowHeight() const;
	
	void setWindowWidth(int w);
	void setWindowHeight(int h);
	
	const std::shared_ptr<TwoDimensionalDisplayController>& getDC() const;
	
	void setCenterX( double x );
	void setCenterY( double y );
	void setCenterZ( double z );
	void setScaleFactor( double scale );
	
	void setCamera( const Camera& cam );
	void setView( const renderingutils::Viewport& view );
    
    int currentCameraIndex() const;

    void finalInit();
	
	void printDDA();
    
    const LiquidInfo& getLiquidInfo();
    
private:
	std::shared_ptr<WetClothCore> m_core;
    std::shared_ptr<TwoDSceneRenderer> m_scene_renderer;
	std::shared_ptr<TwoDimensionalDisplayController> m_display_controller;
	
    TwoDSceneSerializer m_scene_serializer;
};

#endif
