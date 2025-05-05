# CMakeList.txt : CMake project for GasSeger, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (GasSeger "GasSeger.cpp" "Calc.cpp" "GasSeger.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET GasSeger PROPERTY CXX_STANDARD 20)
endif()
#include(src/3rdparty/opennurbs/openNURBS.cmake)
#include(src/on_nurbs/on_nurbs.cmake)
#include(src/3rdparty/opennurbs/zlib.cmake)
#list(APPEND OPENNURBS_INCLUDES ${ZLIB_INCLUDES})
#list(APPEND OPENNURBS_SOURCES ${ZLIB_SOURCES})
#add_definitions(-DNOMINMAX -DNOGDI)

#set(CMAKE_FIND_USE_CMAKE_SYSTEM_PATH OFF)
#set(CMAKE_PREFIX_PATH "E:/repos/src/build")
#cmake_policy(SET CMP0074 NEW)
#set(PCL_ROOT "E:/repos/src/build")
#set(pcl_lib "C:/Program Files/PCL 1.13.0/lib")
find_package(PCL REQUIRED )
#message(STATUS "HERE")
#message(STATUS ${PCL_INCLUDE_DIRS})

#include(src/3rdparty/opennurbs/openNURBS.cmake)
#include(src/on_nurbs/on_nurbs.cmake)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
link_libraries(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
#include_directories(pcl)
include_directories(${VTK_INCLUDE_DIRS})
link_directories(${VTK_LIBRARY_DIRS})
add_definitions(${VTK_DEFINITIONS})
list(REMOVE_ITEM PCL_LIBRARIES pcl_surface)
target_link_libraries (GasSeger PUBLIC ${PCL_LIBRARIES})
set(ENV_PATH "PATH=C:/Program Files/PCL 1.13.0/bin;C:/Program Files/PCL 1.13.0/3rdParty/VTK/bin;C:/Program Files/PCL 1.13.0/3rdParty/flann/bin;")
set_property(TARGET GasSeger PROPERTY
	VS_DEBUGGER_ENVIRONMENT
	${ENV_PATH})
target_link_libraries (GasSeger PUBLIC ${VTK_LIBRARIES})
target_link_libraries (GasSeger PUBLIC OPN)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
# target_link_libraries (GasSeger PUBLIC pcl_common pcl_kdtree pcl_search pcl_features pcl_io pcl_segmentation pcl_sample_consensus)
#message(STATUS "HERE")
#message(STATUS ${PCL_LIBRARIES})
# TODO: Add tests and install targets if needed.
