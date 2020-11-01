# libWetCloth Library

append_files (Headers "h" . ThinShell ThinShell/Forces DER DER/Forces DER/Dependencies pcgsolver)
append_files (Sources "cpp" . ThinShell ThinShell/Forces DER DER/Forces DER/Dependencies pcgsolver)

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

add_library (WetCloth STATIC ${Headers} ${Sources})
target_link_libraries (WetCloth ${LIBWETCLOTH_LIBRARIES})
set_target_properties(WetCloth PROPERTIES
    CXX_STANDARD 14
    CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
