#ifndef __TWO_D_SCENE_XML_PARSER_H__
#define __TWO_D_SCENE_XML_PARSER_H__

#include <Eigen/StdVector>

#include <iostream>
#include <fstream>
#include <limits>

#include "rapidxml.hpp"

#include "TwoDScene.h"

#include "LinearizedImplicitEuler.h"
#include "LevelSetForce.h"

#include "SpringForce.h"
#include "SimpleGravityForce.h"

#include "StringUtilities.h"
#include "RenderingUtilities.h"

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

#include "TwoDSceneRenderer.h"
#include "TwoDimensionalDisplayController.h"
#include "ParticleSimulation.h"

#include "DER/StrandForce.h"
#include "DER/StrandParameters.h"

#include "DistanceFields.h"

#include "ThinShell/ThinShellForce.h"
#include "CohesionForce.h"
#include "JunctionForce.h"

#include "Camera.h"


// REALLY USEFULL TODOs
//   TODO: Improve error messages to display all valid options, etc. Could define an option class that knows its valid options and bounds on values.

// LESS USEFULL TODOs
//   TODO: Write method for computing number of a given property
//   TODO: Add some additional error checking for repeat properties, etc
//   TODO: Abstract out common code
//   TODO: Check for invalid properties

class TwoDSceneXMLParser
{
public:
    
	void loadExecutableSimulation( const std::string& file_name, bool rendering_enabled, std::shared_ptr<ParticleSimulation>& execsim, Camera& cam, scalar& dt, scalar& max_time, scalar& steps_per_sec_cap, renderingutils::Color& bgcolor, std::string& description, std::string& scenetag, bool& cam_inited, const std::string& input_bin );
    
    // TODO: NEED AN EIGEN_ALIGNED_THING_HERE ?
private:
    
	void loadParticleSimulation( bool rendering_enabled, std::shared_ptr<ParticleSimulation>& execsim, scalar& dt, renderingutils::Color& bgcolor, rapidxml::xml_node<>* node, const std::string& input_bin );
    
	void loadLiquidInfo(rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene);
	
    void loadXMLFile( const std::string& filename, std::vector<char>& xmlchars, rapidxml::xml_document<>& doc );
	
    bool loadTextFileIntoString( const std::string& filename, std::string& filecontents );
    
    void loadSimulationType( rapidxml::xml_node<>* node, std::string& simtype );
    
    void loadHairs(rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene, const scalar& dt);
	
	void loadClothes(rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene);

	void loadSpringForces( rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene );
	
	void loadSimpleGravityForces( rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene );
	
    void loadParticles(rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene, int& maxgroup );
    
    void loadSceneTag( rapidxml::xml_node<>* node, std::string& scenetag );
	
    void loadDistanceFields( rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene, int& maxgroup );
	
    void loadStrandParameters( rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene, const scalar& dt );
    
	void loadBucketInfo( rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene );
	
	void loadIntegrator( rapidxml::xml_node<>* node, std::shared_ptr<SceneStepper>& scenestepper, scalar& dt );
    
	void loadScripts( rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene );
	
	bool loadCamera( rapidxml::xml_node<>* node, Camera& camera );
	
    void loadMaxTime( rapidxml::xml_node<>* node, scalar& max_t );
    
    void loadMaxSimFrequency( rapidxml::xml_node<>* node, scalar& max_freq );
    
    void loadViewport( rapidxml::xml_node<> *node, renderingutils::Viewport &view);
    
    void loadBackgroundColor( rapidxml::xml_node<>* node, renderingutils::Color& color );
    
    void loadSceneDescriptionString( rapidxml::xml_node<>* node, std::string& description_string );
};

#endif
