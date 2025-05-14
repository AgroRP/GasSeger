# CMakeList.txt : CMake project for GasSeger, include source and define
# project specific logic here.
#

include(src/3rdparty/opennurbs/openNURBS.cmake)
include(src/on_nurbs/on_nurbs.cmake)
  
include(src/3rdparty/opennurbs/zlib.cmake)
list(APPEND OPENNURBS_INCLUDES ${ZLIB_INCLUDES})
list(APPEND OPENNURBS_SOURCES ${ZLIB_SOURCES})
add_definitions(-DNOMINMAX -DNOGDI)
add_library(OPN STATIC  ${OPENNURBS_SOURCES} ${ON_NURBS_SOURCES})

# Add source to this project's executable.
add_executable (ReConstruct "NurbsFit.cpp")

set_property(TARGET OPN PROPERTY CXX_STANDARD 17)
set_property(TARGET ReConstruct PROPERTY CXX_STANDARD 17)
set(CXXFLAGS "${CXXFLAGS} -march=native")


# TODO:: change this to your PCL lib
set(pcl_lib "C:/Program Files/PCL 1.13.0/lib")
find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
link_libraries(${PCL_LIBRARY_DIRS})
link_directories(libs)
add_definitions(${PCL_DEFINITIONS})
#include_directories(pcl)
include_directories(${VTK_INCLUDE_DIRS})
link_directories(${VTK_LIBRARY_DIRS})
add_definitions(${VTK_DEFINITIONS})
target_link_libraries(OPN pcl_common pcl_search pcl_kdtree pcl_octree ${ON_NURBS_LIBRARIES} RpcRT4.lib)
target_link_libraries (ReConstruct PUBLIC ${pcl_lib}/pcl_gpu_segmentationd.lib ${pcl_lib}/pcl_gpu_octreed.lib ${pcl_lib}/pcl_gpu_containersd.lib)
target_link_libraries (ReConstruct PUBLIC ${VTK_LIBRARIES}  ${ON_NURBS_LIBRARIES} ${PCL_LIBRARIES} OPN)
set(CMAKE_STATIC_LIBRARY_LINK_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_FLAGS} /FORCE:MULTIPLE")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
# target_link_libraries (GasSeger PUBLIC pcl_common pcl_kdtree pcl_search pcl_features pcl_io pcl_segmentation pcl_sample_consensus)
# TODO: Add tests and install targets if needed.
