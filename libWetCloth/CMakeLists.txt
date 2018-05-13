# libWetCloth Executable

append_files (Headers "h" . App Core Core/ThinShell Core/ThinShell/Forces Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)
append_files (Sources "cpp" . App Core Core/ThinShell Core/ThinShell/Forces Core/DER Core/DER/Forces Core/DER/Dependencies Core/pcgsolver)

include_directories (Core)

option (USE_OPENGL "Builds in support for OpenGL rendering" ON)
if (USE_OPENGL)
add_definitions (-DRENDER_ENABLED)
# Locate OpenGL
find_package (OpenGL REQUIRED)
if (OPENGL_FOUND)
  include_directories (${OPENGL_INCLUDE_DIR})
  set (LIBWETCLOTH_LIBRARIES ${LIBWETCLOTH_LIBRARIES} ${OPENGL_LIBRARIES})
else (OPENGL_FOUND)
  message (SEND_ERROR "Unable to locate OpenGL")
endif (OPENGL_FOUND)

# Locate GLUT
find_package (GLUT REQUIRED glut)
if (GLUT_FOUND)
  include_directories (${GLUT_INCLUDE_DIR})
  set (LIBWETCLOTH_LIBRARIES ${LIBWETCLOTH_LIBRARIES} ${GLUT_glut_LIBRARY})
else (GLUT_FOUND)
  message (SEND_ERROR "Unable to locate GLUT")
endif (GLUT_FOUND)

find_package (PNG)
if (PNG_FOUND)
  add_definitions (-DPNGOUT)
  add_definitions (${PNG_DEFINITIONS})
  include_directories (${PNG_INCLUDE_DIR})
  set (LIBWETCLOTH_LIBRARIES ${LIBWETCLOTH_LIBRARIES} ${PNG_LIBRARIES})
endif (PNG_FOUND)

# Locate AntTweakBar
find_package (ANTTWEAKBAR REQUIRED)
if (ANTTWEAKBAR_FOUND)
  include_directories( ${ANT_TWEAK_BAR_INCLUDE_DIR} )
  set (LIBWETCLOTH_LIBRARIES ${LIBWETCLOTH_LIBRARIES} ${ANT_TWEAK_BAR_LIBRARY})  
else (ANTTWEAKBAR_FOUND)
  message (SEND_ERROR "Unable to locate ANTTWEAKBAR")
endif (ANTTWEAKBAR_FOUND)  

endif (USE_OPENGL)

find_package (TBB REQUIRED)
if (TBB_FOUND)
  include_directories (${TBB_INCLUDE_DIRS})
  set (LIBWETCLOTH_LIBRARIES ${LIBWETCLOTH_LIBRARIES} ${TBB_LIBRARIES})
else (TBB_FOUND)
  message (SEND_ERROR "Unable to locate TBB")
endif (TBB_FOUND)

# IGL library is required
find_package (LIBIGL REQUIRED)
if (LIBIGL_FOUND)
  include_directories (${LIBIGL_INCLUDE_DIR})
else (LIBIGL_FOUND)
  message (SEND_ERROR "Unable to locate IGL")
endif (LIBIGL_FOUND)

# RapidXML library is required
find_package (RapidXML REQUIRED)
if (RAPIDXML_FOUND)
  include_directories (${RAPIDXML_INCLUDE_DIR})
else (RAPIDXML_FOUND)
  message (SEND_ERROR "Unable to locate RapidXML")
endif (RAPIDXML_FOUND)

# TCLAP library is required
find_package (TCLAP REQUIRED)
if (TCLAP_FOUND)
  include_directories (${TCLAP_INCLUDE_PATH})
else (TCLAP_FOUND)
  message (SEND_ERROR "Unable to locate TCLAP")
endif (TCLAP_FOUND)

#message(STATUS "Extra libs in libWetCloth: ${LIBWETCLOTH_LIBRARIES}")
#message(STATUS "INSTALL: $CMAKE_INSTALL_PREFIX}")

add_executable (libWetCloth ${Headers} ${Templates} ${Sources})
target_link_libraries (libWetCloth ${LIBWETCLOTH_LIBRARIES})

INSTALL_TARGETS(/bin libWetCloth)
