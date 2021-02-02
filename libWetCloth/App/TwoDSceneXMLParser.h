//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TWO_D_SCENE_XML_PARSER_H
#define TWO_D_SCENE_XML_PARSER_H

#include <Eigen/StdVector>
#include <fstream>
#include <iostream>
#include <limits>

#include "Camera.h"
#include "CohesionForce.h"
#include "DER/StrandForce.h"
#include "ElasticParameters.h"
#include "DistanceFields.h"
#include "JunctionForce.h"
#include "LevelSetForce.h"
#include "LinearizedImplicitEuler.h"
#include "ParticleSimulation.h"
#include "RenderingUtilities.h"
#include "SimpleGravityForce.h"
#include "SpringForce.h"
#include "StringUtilities.h"
#include "ThinShell/ThinShellForce.h"
#include "TwoDScene.h"
#include "TwoDSceneRenderer.h"
#include "TwoDimensionalDisplayController.h"
#include "rapidxml.hpp"

// REALLY USEFULL TODOs
//   TODO: Improve error messages to display all valid options, etc. Could
//   define an option class that knows its valid options and bounds on values.

// LESS USEFULL TODOs
//   TODO: Write method for computing number of a given property
//   TODO: Add some additional error checking for repeat properties, etc
//   TODO: Abstract out common code
//   TODO: Check for invalid properties

class TwoDSceneXMLParser {
 public:
  void loadExecutableSimulation(const std::string& file_name,
                                bool rendering_enabled,
                                std::shared_ptr<ParticleSimulation>& execsim,
                                Camera& cam, scalar& dt, scalar& max_time,
                                scalar& steps_per_sec_cap,
                                renderingutils::Color& bgcolor,
                                std::string& description, std::string& scenetag,
                                bool& cam_inited, const std::string& input_bin);

  // TODO: NEED AN EIGEN_ALIGNED_THING_HERE ?
 private:
  void loadParticleSimulation(bool rendering_enabled,
                              std::shared_ptr<ParticleSimulation>& execsim,
                              scalar& dt, renderingutils::Color& bgcolor,
                              rapidxml::xml_node<>* node,
                              const std::string& input_bin);

  void loadLiquidInfo(rapidxml::xml_node<>* node,
                      const std::shared_ptr<TwoDScene>& twodscene);

  void loadXMLFile(const std::string& filename, std::vector<char>& xmlchars,
                   rapidxml::xml_document<>& doc);

  bool loadTextFileIntoString(const std::string& filename,
                              std::string& filecontents);

  void loadSimulationType(rapidxml::xml_node<>* node, std::string& simtype);

  void loadHairs(rapidxml::xml_node<>* node,
                 const std::shared_ptr<TwoDScene>& twodscene, const scalar& dt);

  void loadHairPose(rapidxml::xml_node<>* node,
                    const std::shared_ptr<TwoDScene>& twodscene);

  void loadClothes(rapidxml::xml_node<>* node,
                   const std::shared_ptr<TwoDScene>& twodscene);

  void loadSpringForces(rapidxml::xml_node<>* node,
                        const std::shared_ptr<TwoDScene>& twodscene);

  void loadSimpleGravityForces(rapidxml::xml_node<>* node,
                               const std::shared_ptr<TwoDScene>& twodscene);

  void loadParticles(rapidxml::xml_node<>* node,
                     const std::shared_ptr<TwoDScene>& twodscene,
                     int& maxgroup);

  void loadSceneTag(rapidxml::xml_node<>* node, std::string& scenetag);

  void loadDistanceFields(rapidxml::xml_node<>* node,
                          const std::shared_ptr<TwoDScene>& twodscene,
                          int& maxgroup);

  void loadElasticParameters(rapidxml::xml_node<>* node,
                            const std::shared_ptr<TwoDScene>& twodscene,
                            const scalar& dt);

  void loadBucketInfo(rapidxml::xml_node<>* node,
                      const std::shared_ptr<TwoDScene>& twodscene);

  void loadIntegrator(rapidxml::xml_node<>* node,
                      std::shared_ptr<SceneStepper>& scenestepper, scalar& dt);

  void loadScripts(rapidxml::xml_node<>* node,
                   const std::shared_ptr<TwoDScene>& twodscene);

  bool loadCamera(rapidxml::xml_node<>* node, Camera& camera);

  void loadMaxTime(rapidxml::xml_node<>* node, scalar& max_t);

  void loadMaxSimFrequency(rapidxml::xml_node<>* node, scalar& max_freq);

  void loadViewport(rapidxml::xml_node<>* node, renderingutils::Viewport& view);

  void loadBackgroundColor(rapidxml::xml_node<>* node,
                           renderingutils::Color& color);

  void loadSceneDescriptionString(rapidxml::xml_node<>* node,
                                  std::string& description_string);
};

#endif
