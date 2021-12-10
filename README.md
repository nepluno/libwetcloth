[![libWetCloth](http://www.cs.columbia.edu/cg/raymond/tighten_the_towel.jpg)](http://www.cs.columbia.edu/cg/wetcloth/)

libWetCloth
================
libWetCloth is an open source project for the physical simulation of liquid and wet cloth or yarns. It has been compiled and tested on Mac OS X (with both Intel and Apple M1 chips), Ubuntu Linux, Windows, and licensed under the Mozilla Public License v. 2.0.

We would like to hear from you if you appreciate this work.

It is the original implementation of paper A Multi-Scale Model for Simulating Liquid-Fabric Interactions ( refer to our project page for more details: http://libWetCloth.info, or http://www.cs.columbia.edu/cg/wetcloth/ ). This code base contains the following parts:

 - A liquid simulator implementing the affine-particle-in-cell method.
 - A cloth simulator implementing the elastic thin shell model.
 - A yarn simulator implementing the discrete elastic rod (DER) / viscous thread model.
 - A cloth/yarn collision handler based on anisotropic elastoplasticity, discretized with augmented, moving least sqaures material point method (AMLS-MPM). 
 - A two-way coupling method based on mixture theory, between the cloth, yarn and liquid, handling dragging, buoyancy, capturing and dripping effect.

Directories
--------------------
The directories of this code repository are structured as following:
```
libwetcloth
│
└───assets: files to define scenarios
│   │
│   └───buoyancy_tests: scenarios to test buoyancy effects
│   │
│   └───drag_tests: scenarios to test different drag models
│   │
│   └───general_examples: scenarios used in the paper as demos
│   │
│   └───parameter_tests: scenarios to test various liquid materials 
│   │                    and yarn settings
│   │
│   └───ring_tests: the ring test scenario used in the paper
│   │
│   └───unit_tests: simple scenarios to test dynamics and collisions 
│                   (including dry yarn/cloth)
│
└───cmake: CMake files to find packages
│
└───houdini: Houdini scenes to produce renderings used in the paper
│
└───include: thirdparty libraries (Eigen, libIGL, etc.)
│
└───libWetCloth
    │
    └───App: user-interface to interact with simulation
    │
    └───Core: actual simulation code
        │
        └───DER: code for the discrete elastic rods (for yarns/hairs)
        │   │
        │   └───Dependencies: computing the geometric information
        │   │                 used in the discrete elastic rods
        │   │                 model, including the reference/material
        │   │                 frames, reference/material twists, and
        │   │                 the material curvatures.
        │   │
        │   └───Forces: computing the stretching, bending, and twisting
        │               forces and their Jacobians.
        │
        └───ThinShell: code for the discrete elastic thin shell
        │   │          (for triangular clothes) 
        │   │
        │   └───Forces: computing the stretching and bending forces, as
        │               well as their Jacobians.
        │
        └───PCGSolver: PCG solver with the MIC preconditioner
```

Dependencies
--------------------
libWetCloth depends on following libraries (some of them have been included in the code base for all platforms, marked with an asterisk):

- Eigen* (http://eigen.tuxfamily.org/)
- RapidXML* (http://rapidxml.sourceforge.net)
- tclap* (http://tclap.sourceforge.net)
- libIGL* (https://github.com/libigl/libigl)
- AntTweakBar (http://anttweakbar.sourceforge.net/doc/)
- Intel TBB (https://www.threadingbuildingblocks.org)
- FreeGLUT (http://freeglut.sourceforge.net/)
- libPNG (https://libpng.sourceforge.io/)
- zlib (https://www.zlib.net/)

On Mac OS X or Linux-based systems, most of the dependencies are either included, or can be easily installed with Homebrew (https://brew.sh) or the APT package handling utility. For Intel TBB, you may download and install from the link provided above, or from Intel website (https://software.intel.com/en-us/parallel-studio-xe/choose-download). 

On Windows you may need manually download and compile some of them (e.g. AntTweakBar, TBB, libPNG). For the ease of compilation, we provide a package containing all the headers and pre-compiled thirdparty libraries (except for those that have been included for all platforms). Please refer to the `On Windows` section below for more details.

You may also compile a version that directly run simulation in the console, by turning off the `USE_OPENGL` switch in the CMake settings. When this switch is turned off, dependencies such as **GLUT, libPNG, zlib, and AntTweakBar** are no longer necessary to compile libWetCloth.

For more details, please refer to the compilation section below.

Compilation
-----------------
libWetCloth has been tested with Clang (under Mac OS X), GCC 4.8+ (under Linux), and Microsoft Visual Studio (under Windows 10).

To compile libWetCloth, you'll need CMake on Mac OS X or Linux, or CMake-GUI (https://cmake.org) on Windows.

On Mac OS X or Linux:

1. make a `<build>` directory, say, just *build*, with *mkdir build*, enter the `<build>` directory, type *cmake ..*
2. Optionally you can adjust the options with *ccmake ..* In some cases there can be some packages that cmake cannot find. You need to manually specify their paths through ccmake then.
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

On Windows:

You may download the pre-compiled package (https://www.cs.columbia.edu/cg/raymond/libwetcloth_thirdparty_win64.zip) that contains the x64 pre-compiled binaries for AntTweakBar, Intel TBB, GLUT, libPNG and zlib, unpack it somewhere (e.g., `<libwetcloth directory>/thirdparty`), and then specify the missing directories to the path containing the headers or the compiled libraries. If you compile it manually, please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86, and 64-bit libraries for x64).

If our provided third-party libraries are used, the following CMake variables should be set before configuration (replace the `<libwetcloth_thirdparty_win64>` with the actual unpacked directory path. If `<libwetcloth directory>/thirdparty` is used for `<libwetcloth_thirdparty_win64>`, the dependencies should be found by CMake automatically):
- `ANT_TWEAK_BAR_INCLUDE_DIR`: `<libwetcloth_thirdparty_win64>/include/AntTweakBar`
- `ANT_TWEAK_BAR_LIBRARY`: `<libwetcloth_thirdparty_win64>/lib/AntTweakBar64.lib`
- `FREEGLUT_LIBRARY`: `<libwetcloth_thirdparty_win64>/lib/freeglut.lib`
- `PNG_LIBRARY_DEBUG`: `<libwetcloth_thirdparty_win64>/lib/libpng16_staticd.lib`
- `PNG_LIBRARY_RELEASE`: `<libwetcloth_thirdparty_win64>/lib/libpng16_static.lib`
- `PNG_PNG_INCLUDE_DIR`: `<libwetcloth_thirdparty_win64>/include/libpng`
- `TBB_INCLUDE_DIRS`: `<libwetcloth_thirdparty_win64>/include`
- `TBB_tbb_LIBRARY_DEBUG`: `<libwetcloth_thirdparty_win64>/lib/tbb_debug.lib`
- `TBB_tbb_LIBRARY_RELEASE`: `<libwetcloth_thirdparty_win64>/lib/tbb.lib`
- `ZLIB_INLCUDE_DIR`: `<libwetcloth_thirdparty_win64>/include/zlib`
- `ZLIB_LIBRARY_DEBUG`: `<libwetcloth_thirdparty_win64>/lib/zlibstaticd.lib`
- `ZLIB_LIBRARY_RELEASE`: `<libwetcloth_thirdparty_win64>/lib/zlibstatic.lib`

1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*, choose your installed version of the Microsoft Visual Studio.
2. after configuration you may find several libraries not found (with notifications of errors), check the *Advanced* box and *specify those missing header path and libraries manually*. For example, if Eigen is missing, then please specify the EIGEN3_INCLUDE_DIR to the path of directory we provided. 
3. click generate after fixing all missing variables to generate your Visual Studio solution.
4. open the Visual Studio solution and compile the code.
5. before running the demo, all the compiled dynamic linking libraries (DLLs) for your dependencies should be accessible from the executable. The DLLs should be either in the same directory with the executable, or your PATH environment variable that can be changed in system settings. You may copy the DLLs (we include them in `<libwetcloth_thirdparty_win64>/bin` as the pre-compiled version) to the path of the compiled executable (`<build>/libWetCloth/App/Release` or `<build>/libWetCloth/App/Debug`) or you may simply copy them into your System32 (x64) or SysWOW64 (x86) directories.

Run the Demo
--------------------
To run the demo of libWetCloth, you may simply use the command line argument *-s [scene_file]* to specify the scene to be loaded. Besides of the WetCloth library, we also have a user interface named WetClothApp where you may watch and tune the simulation interactively. For example, you may type
```
./WetClothApp -s ../assets/unit_tests/simple_cloth.xml
```
under the `<build>/libWetCloth/App` directory to run the simulation of the scene containing a water ball splashes on a small cloth. 

All the parameters can be modified offline in the scene description XML files. Some can be changed online in the user interface provided by the demo program.

USAGE: 
```
   ./WetClothApp -s <string> [-i <string>] [-o <integer>] [-g <integer>] [-d <boolean>] [-p <boolean>] [--] [--version] [-h]

Where: 
   -s <string>,  --scene <string>
     (required)  Simulation to run; an xml scene file

   -i <string>,  --inputfile <string>
     Binary file to load simulation pos from

   -o <integer>,  --outputfile <integer>
     Binary file to save simulation state to

   -g <integer>,  --generate <integer>
     Generate PNG if 1, not if 0

   -d <boolean>,  --display <boolean>
     Run the simulation with display enabled if 1, without if 0

   -p <boolean>,  --paused <boolean>
     Begin the simulation paused if 1, running if 0

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

Surface Reconstruction and Rendering with Houdini
--------------------------------------------------------
The Houdini projects are also provided in the "houdini" folder, which are used for surface reconstruction and rendering purposes. Our simulator can generate data that can be read back by the Python script in our Houdini projects.

You may run the demo with the "-o" option, *under the folder containing Houdini projects* (by default, [project source]/houdini/). For example, you may type
```
../build/libWetCloth/App/WetClothApp -s ../assets/parameter_tests/pore_test_mid.xml -o 20
```
under the folder *containing Houdini projects* to generate data per 20 time steps (for this example we use 0.0002s for the simulation time step, and 0.004s for the rendering time step. Hence 20 time steps is taken for the data generation). The simulation code will create a folder with the name of the scene file under this folder ([project source]/houdini/pore_test_mid for this example), which can be read back by the Houdini project with the same name.

After some data generated, you may open the corresponding Houdini project to watch and operate on them. We use the nodes with suffix "bake" to indicate the usage of baking. For example in the "fluids_bake" node, the fluid particles will be read and used to reconstruct a polygonal liquid surface, which will then be stored as Houdini geometry files.

The baked geometries are then read back by the nodes with prefix "geo", which can be used for rendering. 

You may check the Mantra nodes in the output classifier for the rendering settings. You can also connect these nodes to HQueue nodes for distributed rendering on a render farm.

Houdini 16.5 is used for our case. You may need an equivalent or higher version to open the project files. For more information and tutorials of Houdini, please visit the SideFX website (https://www.sidefx.com/).

Contact
-----------
Please contact the author (yf2320@columbia.edu) for questions and bug report, or consulting for the usage of this code base.

BibTex Citation
----------------------
@article{Fei:2018:MMS:3197517.3201392  
 author = {Fei, Yun (Raymond) and Batty, Christopher and Grinspun, Eitan and Zheng, Changxi},  
 title = {A Multi-scale Model for Simulating Liquid-fabric Interactions},  
 journal = {ACM Trans. Graph.},  
 issue_date = {Aug 2018},  
 volume = {37},  
 number = {4},  
 month = aug,  
 year = {2018},  
 pages = {51:1--51:16},  
 articleno = {51},  
 numpages = {16},  
 url = {http://doi.acm.org/10.1145/3197517.3201392},  
 doi = {10.1145/3197517.3201392},  
 acmid = {3201392},  
 publisher = {ACM},  
 address = {New York, NY, USA},  
 keywords = {fluid dynamics, particle-in-cell, mixture theory, two-way coupling, wet cloth},  
}  