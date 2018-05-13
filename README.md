libWetCloth
================
libWetCloth is an open source project for the physical simulation of liquid and wet cloth or yarns. It is cross-platform (Mac OS X, Linux, Windows, and more), and licensed under Clear BSD License for academic and non-commercial use (other licenses may be obtained by contacting the faculty of the Columbia Computer Graphics Group or a Columbia University licensing officer). 

We would like to hear from you if you appreciate this work.

It is the original implementation of paper A Multi-Scale Model for Simulating Liquid-Fabric Interactions ( refer to our project page for more details: http://libWetCloth.info, or http://www.cs.columbia.edu/cg/wetcloth/ ). This code base contains the following parts:

 - A liquid simulator implementing the affine-particle-in-cell method.
 - A cloth simulator implementing the elastic thin shell model.
 - A yarn simulator implementing the elastic viscous thread model.
 - An MPM-based collision handler (for more details, see Anisotropic Elastoplasticity for Cloth, Knit and Hair Frictional Contact, http://www.seas.upenn.edu/~cffjiang/research/cloth/paper.pdf)
 - A mixture theory based two-way coupling method between the cloth, yarn and liquid, including dragging, buoyancy, capturing and dripping effect.

Dependencies
--------------------
libWetCloth depends on following libraries (some of them have been included in the code base, marked with asterisk):

- Eigen* (http://eigen.tuxfamily.org/)
- RapidXML* (http://rapidxml.sourceforge.net)
- tclap* (http://tclap.sourceforge.net)
- libIGL* (https://github.com/libigl/libigl)
- AntTweakBar (http://anttweakbar.sourceforge.net/doc/)
- Intel TBB (https://www.threadingbuildingblocks.org)
- GLUT (https://www.opengl.org/resources/libraries/glut/)
- libPNG (https://libpng.sourceforge.io/)

On Mac OS X or Linux-based systems, most of the dependencies are either included, or can be easily installed with Homebrew (https://brew.sh) or the APT package handling utility. For Intel TBB, you may download and install from the link provided above, or from Intel website (https://software.intel.com/en-us/parallel-studio-xe/choose-download).

On Windows you may need manually download and install some of them.

Compilation
-----------------
libWetCloth has been tested with Clang (under Mac OS X), GCC 4.8+ (under Linux), and Microsoft Visual Studio (under Windows 10).

To compile libWetCloth, you'll need CMake on Mac OS X or Linux, or CMake-GUI (https://cmake.org) on Windows.

On Mac OS X or Linux:
1. make a directory, say, *build*, with *mkdir build*, enter the *build* directory, type *cmake ..*
2. Optionally you can adjust the options with *ccmake .*
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

On Windows:
1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*, choose your installed version of the Microsoft Visual Studio.
2. after configuration you may find several libraries not found, check the *Advanced* box and locate those missing libraries manually. Please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86, and 64-bit libraries for x64).
3. click generate after fixing all missing variables to generate your Visual Studio solution.
4. open the solution and compile the code.

Run the Demo
--------------------
To run the demo of libWetCloth, you may simply use the command line argument *-s [scene_file]* to specify the scene to be loaded. For example, you may type

./libWetCloth -s assets/general_examples/splash_cloth_small.xml

to run the simulation of the scene containing a water ball splashes on a small cloth. 

All the parameters can be modified offline in the scene description XML files. Some can be changed online in the user interface provided by the demo program.

USAGE: 

   ./libWetCloth -s <string> [-i <string>] [-o <integer>] [-g <integer>] [-d <boolean>] [-p <boolean>] [--] [--version] [-h]


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