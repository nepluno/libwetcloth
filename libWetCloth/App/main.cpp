//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sys/stat.h>
#include <tclap/CmdLine.h>

#include <Eigen/StdVector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>

#ifdef WIN32
#include <direct.h>

#ifndef GL_MULTISAMPLE_ARB
#define GL_MULTISAMPLE_ARB 0x809D
#endif

#ifndef GL_MULTISAMPLE_FILTER_HINT_NV
#define GL_MULTISAMPLE_FILTER_HINT_NV 0x8534
#endif
#endif

#include "Camera.h"
#include "Force.h"
#include "MathDefs.h"
#include "StringUtilities.h"
#include "TimingUtilities.h"
#include "TwoDScene.h"
#include "TwoDSceneSerializer.h"
#include "TwoDSceneXMLParser.h"

#ifdef RENDER_ENABLED
#include <AntTweakBar.h>

#include "RenderingUtilities.h"
#include "TwoDSceneRenderer.h"
#include "TwoDimensionalDisplayController.h"
#include "YImage.h"
#endif

#include "ParticleSimulation.h"

///////////////////////////////////////////////////////////////////////////////
// Contains the actual simulation, renderer, parser, and serializer
std::shared_ptr<ParticleSimulation> g_executable_simulation = NULL;

///////////////////////////////////////////////////////////////////////////////
// Rendering State
bool g_rendering_enabled = true;
double g_sec_per_frame;
double g_last_time = timingutils::seconds();

renderingutils::Color g_bgcolor(1.0, 1.0, 1.0);

///////////////////////////////////////////////////////////////////////////////
// SVG Rendering State

///////////////////////////////////////////////////////////////////////////////
// Parser state
std::string g_xml_scene_file;
std::string g_description;
std::string g_scene_tag = "";

///////////////////////////////////////////////////////////////////////////////
// Scene input/output/comparison state
int g_save_to_binary = 0;
std::string g_binary_file_name;
std::ofstream g_binary_output;
std::string g_short_file_name;

///////////////////////////////////////////////////////////////////////////////
// Simulation state
int g_dump_png = 0;
bool g_paused = true;
scalar g_dt = 0.0;
int g_num_steps = 0;
int g_current_step = 0;

///////////////////////////////////////////////////////////////////////////////
// Simulation functions

void miscOutputCallback();
void dumpPNG(const std::string& filename);

void stepSystem() {
  g_executable_simulation->stepSystem(g_dt);
  g_current_step++;

  // Execute the user-customized output callback
  miscOutputCallback();

  // Determine if the simulation is complete
  if (g_current_step >= g_num_steps) {
    std::cout << "Complete Simulation! Enter time to continue (exit with 0): "
              << std::endl;
    double new_time = 0.0;
    std::cin >> new_time;

    g_num_steps += ceil(new_time / g_dt);
    if (new_time < g_dt) g_paused = true;
  }
}

void syncScene() {}

void headlessSimLoop() {
  while (true) {
    while (g_current_step <= g_num_steps) {
      stepSystem();
    }

    std::cout << "Complete Simulation! Enter time to continue (exit with 0): "
              << std::endl;
    double new_time = 0.0;
    std::cin >> new_time;
    if (new_time < g_dt) break;

    g_num_steps += ceil(new_time / g_dt);
  }
}

#ifdef RENDER_ENABLED
void dumpPNGsubprog(YImage* image, char* fnstr) {
  image->flip();
  image->save(fnstr);
  delete image;
  delete fnstr;
}

///////////////////////////////////////////////////////////////////////////////
// Rendering and UI functions

void dumpPNG(const std::string& filename) {
  YImage* image = new YImage;
  char* fnstr = new char[filename.length() + 1];
  strcpy(fnstr, filename.data());

  image->resize(g_executable_simulation->getWindowWidth(),
                g_executable_simulation->getWindowHeight());

  glFinish();

  glPixelStorei(GL_PACK_ALIGNMENT, 4);
  glPixelStorei(GL_PACK_ROW_LENGTH, 0);
  glPixelStorei(GL_PACK_SKIP_ROWS, 0);
  glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
  glReadBuffer(GL_BACK);

  glFinish();
  glReadPixels(0, 0, g_executable_simulation->getWindowWidth(),
               g_executable_simulation->getWindowHeight(), GL_RGBA,
               GL_UNSIGNED_BYTE, image->data());

  std::thread t(std::bind(dumpPNGsubprog, image, fnstr));

  t.detach();
}

void reshape(int w, int h) {
  TwWindowSize(w, h);

  g_executable_simulation->reshape(w, h);

  assert(renderingutils::checkGLErrors());
}

// TODO: Move these functions to scene renderer?
void setOrthographicProjection() {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  gluOrtho2D(0, g_executable_simulation->getWindowWidth(), 0,
             g_executable_simulation->getWindowHeight());

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  assert(renderingutils::checkGLErrors());
}

void renderBitmapString(float x, float y, float z, void* font, std::string s) {
  glRasterPos3f(x, y, z);
  for (std::string::iterator i = s.begin(); i != s.end(); ++i) {
    char c = *i;
    glutBitmapCharacter(font, c);
  }

  assert(renderingutils::checkGLErrors());
}

void drawHUD() {
  setOrthographicProjection();
  glColor3f(1.0 - g_bgcolor.r, 1.0 - g_bgcolor.g, 1.0 - g_bgcolor.b);
  renderBitmapString(4, g_executable_simulation->getWindowHeight() - 20, 0.0,
                     GLUT_BITMAP_HELVETICA_18,
                     stringutils::convertToString(g_current_step * g_dt));
  renderBitmapString(4, g_executable_simulation->getWindowHeight() - 50, 0.0,
                     GLUT_BITMAP_HELVETICA_18,
                     std::string("Camera [") +
                         stringutils::convertToString(
                             g_executable_simulation->currentCameraIndex()) +
                         std::string("]"));
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  TwDraw();

  assert(renderingutils::checkGLErrors());
}

void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  assert(renderingutils::checkGLErrors());

  assert(g_executable_simulation != NULL);
  g_executable_simulation->renderSceneOpenGL(g_dt);

  drawHUD();

  glutSwapBuffers();

  assert(renderingutils::checkGLErrors());
}

void keyboard(unsigned char key, int x, int y) {
  if (TwEventKeyboardGLUT(key, x, y)) {
    glutPostRedisplay();
    return;
  }

  g_executable_simulation->keyboard(key, x, y);

  if (key == 27 || key == 'q') {
    TwTerminate();

    exit(0);
  } else if (key == 's' || key == 'S') {
    stepSystem();
    glutPostRedisplay();
  } else if (key == ' ') {
    g_paused = !g_paused;
  } else if (key == 'c' || key == 'C') {
    g_executable_simulation->centerCamera();

    glutPostRedisplay();
  } else if (key == 'a' || key == 'A') {
    g_executable_simulation->printDDA();
  }

  assert(renderingutils::checkGLErrors());
}

// Proccess 'special' keys
void special(int key, int x, int y) {
  if (TwEventSpecialGLUT(key, x, y)) {
    glutPostRedisplay();
    return;
  }

  g_executable_simulation->special(key, x, y);

  assert(renderingutils::checkGLErrors());
}

void mouse(int button, int state, int x, int y) {
  if (TwEventMouseButtonGLUT(button, state, x, y)) {
    glutPostRedisplay();
    return;
  }

  g_executable_simulation->mouse(button, state, x, y);

  assert(renderingutils::checkGLErrors());
}

void motion(int x, int y) {
  if (TwEventMouseMotionGLUT(x, y)) {
    glutPostRedisplay();
    return;
  }

  g_executable_simulation->motion(x, y);

  assert(renderingutils::checkGLErrors());
}

void idle() {
  // std::cout << "g_last_time: " << g_last_time << std::endl;
  // Trigger the next timestep
  double current_time = timingutils::seconds();
  // std::cout << "current_time: " << current_time << std::endl;
  // std::cout << "g_sec_per_frame: " << g_sec_per_frame << std::endl;
  if (!g_paused && current_time - g_last_time >= g_sec_per_frame) {
    g_last_time = current_time;
    stepSystem();
    glutPostRedisplay();
  }

#ifdef __APPLE__
  static int mojave_counter = 0;
  if (mojave_counter < 2)
    glutReshapeWindow(
        g_executable_simulation->getWindowWidth() - 1 + mojave_counter++,
        g_executable_simulation->getWindowHeight());
#endif

  assert(renderingutils::checkGLErrors());
}

void initializeOpenGLandGLUT(int argc, char** argv) {
  // Initialize GLUT
  glutInit(&argc, argv);
  // glutInitDisplayString("rgba stencil=8 depth=24 samples=4 hidpi");
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE | GLUT_DEPTH);

  int scrwidth = glutGet(GLUT_SCREEN_WIDTH);
  int scrheight = glutGet(GLUT_SCREEN_HEIGHT);

  int factor = std::min(scrwidth / 320, scrheight / 180);
  scrwidth = factor * 320;
  scrheight = factor * 180;
  /*#ifdef __APPLE__
   scrwidth *= 2;
   scrheight *= 2;
   #endif*/
  g_executable_simulation->setWindowWidth(scrwidth);
  g_executable_simulation->setWindowHeight(scrheight);

  glutInitWindowSize(g_executable_simulation->getWindowWidth(),
                     g_executable_simulation->getWindowHeight());
  glutCreateWindow(
      "libWetCloth: A Multi-Scale Model for Liquid-Fabric Interactions");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutIdleFunc(idle);

  glEnable(GL_MULTISAMPLE_ARB);
  glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);

  TwInit(TW_OPENGL, NULL);
  TwGLUTModifiersFunc(glutGetModifiers);

  // Initialize OpenGL
  reshape(scrwidth, scrheight);
  glClearColor(g_bgcolor.r, g_bgcolor.g, g_bgcolor.b, 1.0);

  glEnable(GL_BLEND);
  glEnable(GL_POLYGON_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

  GLubyte halftone[] = {
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA,
      0x55, 0x55, 0x55, 0x55, 0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55,
      0xAA, 0xAA, 0xAA, 0xAA, 0x55, 0x55, 0x55, 0x55};

  glPolygonStipple(halftone);

  g_executable_simulation->initializeOpenGLRenderer();

  assert(renderingutils::checkGLErrors());
}
#endif

///////////////////////////////////////////////////////////////////////////////
// Parser functions

void loadScene(const std::string& file_name) {
  // Maximum time in the simulation to run for. This has nothing to do with run
  // time, cpu time, etc. This is time in the 'virtual world'.
  scalar max_time;
  // Maximum frequency, in wall clock time, to execute the simulation for. This
  // serves as a cap for simulations that run too fast to see a solution.
  scalar steps_per_sec_cap = 100.0;
  // Contains the center and 'scale factor' of the view
  renderingutils::Viewport view;

  // Load the simulation and pieces of rendring and UI state
  assert(g_executable_simulation == NULL);
  TwoDSceneXMLParser xml_scene_parser;
  bool cam_init = false;
  Camera cam;

  xml_scene_parser.loadExecutableSimulation(
      file_name, g_rendering_enabled, g_executable_simulation, cam, g_dt,
      max_time, steps_per_sec_cap, g_bgcolor, g_description, g_scene_tag,
      cam_init, g_binary_file_name);
  assert(g_executable_simulation != NULL);

  // If the user did not request a custom viewport, try to compute a reasonable
  // default.
  if (!cam_init) {
    g_executable_simulation->centerCamera(false);
  }
  // Otherwise set the viewport per the user's request.
  else {
    g_executable_simulation->setCamera(cam);
  }

  g_executable_simulation->finalInit();

  // To cap the framerate, compute the minimum time a single timestep should
  // take
  g_sec_per_frame = 1.0 / steps_per_sec_cap;
  // Integer number of timesteps to take
  g_num_steps = ceil(max_time / g_dt);
  // We begin at the 0th timestep
  g_current_step = 0;
}

void parseCommandLine(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd(
        "A Multi-Scale Model for Simulating Liquid-Fabric Interactions");

    // XML scene file to load
    TCLAP::ValueArg<std::string> scene("s", "scene",
                                       "Simulation to run; an xml scene file",
                                       true, "", "string", cmd);

    // Begin the scene paused or running
    TCLAP::ValueArg<bool> paused(
        "p", "paused", "Begin the simulation paused if 1, running if 0", false,
        true, "boolean", cmd);

    // Run the simulation with rendering enabled or disabled
    TCLAP::ValueArg<bool> display(
        "d", "display",
        "Run the simulation with display enabled if 1, without if 0", false,
        true, "boolean", cmd);

    // Begin the scene paused or running
    TCLAP::ValueArg<int> dumppng("g", "generate", "Generate PNG if 1, not if 0",
                                 false, 0, "integer", cmd);

    // These cannot be set at the same time
    // File to save output to
    TCLAP::ValueArg<int> output("o", "outputfile",
                                "Binary file to save simulation state to",
                                false, 0, "integer", cmd);
    // File to load for comparisons
    TCLAP::ValueArg<std::string> input(
        "i", "inputfile", "Binary file to load simulation pos from", false, "",
        "string", cmd);

    cmd.parse(argc, argv);

    assert(scene.isSet());
    g_xml_scene_file = scene.getValue();
    g_paused = paused.getValue();
    g_rendering_enabled = display.getValue();
    g_dump_png = dumppng.getValue();
    g_save_to_binary = output.getValue();
    g_binary_file_name = input.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.what() << std::endl;
    exit(1);
  }
}

void cleanupAtExit() {}

std::ostream& main_header(std::ostream& stream) {
  stream << outputmod::startgreen
         << ".__  ._____.   __      __        __  _________ .__          __  "
            ".__     \n"
         << "|  | |__\\_ |__/  \\    /  \\ _____/  |_\\_   ___ \\|  |   _____/ "
            " |_|  |__  \n"
         << "|  | |  || __ \\   \\/\\/   _/ __ \\   __/    \\  \\/|  |  /  _ "
            "\\   __|  |  \\ \n"
         << "|  |_|  || \\_\\ \\        /\\  ___/|  | \\     \\___|  |_(  <_> "
            "|  | |   Y  \\\n"
         << "|____|__||___  /\\__/\\  /  \\___  |__|  \\______  "
            "|____/\\____/|__| |___|  /\n"
         << "             \\/      \\/       \\/             \\/               "
            "       \\/ "
         << outputmod::endgreen << std::endl;

  return stream;
}

std::ofstream g_debugoutput;

void miscOutputCallback() {
  // If the user wants to save output to a binary
  if (g_save_to_binary && !(g_current_step % g_save_to_binary) &&
      g_current_step <= g_num_steps) {
    std::stringstream oss_cloth;
    oss_cloth << g_short_file_name << "/cloth" << std::setw(5)
              << std::setfill('0') << (g_current_step / g_save_to_binary)
              << ".obj";

    std::stringstream oss_hairs;
    oss_hairs << g_short_file_name << "/hair" << std::setw(5)
              << std::setfill('0') << (g_current_step / g_save_to_binary)
              << ".obj";

    std::stringstream oss_fluid;
    oss_fluid << g_short_file_name << "/fluid" << std::setw(5)
              << std::setfill('0') << (g_current_step / g_save_to_binary)
              << ".obj";

    std::stringstream oss_inbd;
    oss_inbd << g_short_file_name << "/internal" << std::setw(5)
             << std::setfill('0') << (g_current_step / g_save_to_binary)
             << ".obj";

    std::stringstream oss_exbd;
    oss_exbd << g_short_file_name << "/external" << std::setw(5)
             << std::setfill('0') << (g_current_step / g_save_to_binary)
             << ".obj";

    std::stringstream oss_spring;
    oss_spring << g_short_file_name << "/spring" << std::setw(5)
               << std::setfill('0') << (g_current_step / g_save_to_binary)
               << ".obj";

    g_executable_simulation->serializeScene(oss_cloth.str(), oss_hairs.str(),
                                            oss_fluid.str(), oss_inbd.str(),
                                            oss_exbd.str(), oss_spring.str());
  }

  // Update the state of the renderers
#ifdef RENDER_ENABLED
  if (g_rendering_enabled) g_executable_simulation->updateOpenGLRendererState();

  // If comparing simulations, load comparison scene's equivalent step

  // If the user wants to generate a PNG movie
  if (g_dump_png && !(g_current_step % g_dump_png)) {
    std::stringstream oss;
    oss << g_short_file_name << "/frame" << std::setw(5) << std::setfill('0')
        << (g_current_step / g_dump_png) << ".png";
    dumpPNG(oss.str());
    oss.clear();
  }
#endif
}

int main(int argc, char** argv) {
  Eigen::initParallel();
  Eigen::setNbThreads(std::thread::hardware_concurrency());

  srand(0x0108170F);

  // Parse command line arguments
  parseCommandLine(argc, argv);

  std::vector<std::string> pathes;

  stringutils::split(g_xml_scene_file, '/', pathes);

  std::vector<std::string> path_first;

  stringutils::split(pathes[pathes.size() - 1], '.', path_first);

  g_short_file_name = path_first[0];

#ifdef WIN32
  _mkdir(g_short_file_name.c_str());
#else
  mkdir(g_short_file_name.c_str(), 0777);
#endif

  // Function to cleanup at progarm exit
  atexit(cleanupAtExit);

  // Load the user-specified scene
  loadScene(g_xml_scene_file);

  // If requested, open the input file for the scene to benchmark
#ifdef RENDER_ENABLED
  // Initialization for OpenGL and GLUT
  if (g_rendering_enabled) initializeOpenGLandGLUT(argc, argv);
#endif
  // Print a header
  std::cout << main_header << std::endl;

#ifdef CMAKE_BUILD_TYPE
  std::cout << outputmod::startblue << "Build type: " << outputmod::endblue
            << CMAKE_BUILD_TYPE << std::endl;
#endif
#ifdef EIGEN_VECTORIZE
  std::cout << outputmod::startblue << "Vectorization: " << outputmod::endblue
            << "Enabled" << std::endl;
#else
  std::cout << outputmod::startblue << "Vectorization: " << outputmod::endblue
            << "Disabled" << std::endl;
#endif

  std::cout << outputmod::startblue << "Scene: " << outputmod::endblue
            << g_xml_scene_file << std::endl;
  std::cout << outputmod::startblue << "Integrator: " << outputmod::endblue
            << g_executable_simulation->getSolverName() << std::endl;

  std::cout << outputmod::startblue
            << "Global Parameters: " << outputmod::endblue << std::endl
            << g_executable_simulation->getLiquidInfo() << std::endl;
#ifdef RENDER_ENABLED
  if (g_rendering_enabled)
    glutMainLoop();
  else
    headlessSimLoop();
#else
  headlessSimLoop();
#endif
  return 0;
}
