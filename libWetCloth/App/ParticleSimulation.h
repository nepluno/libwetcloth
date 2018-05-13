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
