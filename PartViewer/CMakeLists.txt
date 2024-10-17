# CMakeList.txt : CMake project for GasSeger, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (PartsViewer "PartsViewer.cpp" "PartsViewer.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET PartsViewer PROPERTY CXX_STANDARD 20)
endif()


set(pcl_lib "C:/Program Files/PCL 1.13.0/lib")
find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
link_libraries(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
#include_directories(pcl)
include_directories(${VTK_INCLUDE_DIRS})
link_directories(${VTK_LIBRARY_DIRS})
add_definitions(${VTK_DEFINITIONS})
target_link_libraries (PartsViewer PUBLIC ${PCL_LIBRARIES})
target_link_libraries (PartsViewer PUBLIC ${pcl_lib}/pcl_gpu_segmentationd.lib)
target_link_libraries (PartsViewer PUBLIC ${pcl_lib}/pcl_gpu_octreed.lib)
target_link_libraries (PartsViewer PUBLIC ${pcl_lib}/pcl_gpu_containersd.lib)
target_link_libraries (PartsViewer PUBLIC ${VTK_LIBRARIES})
# target_link_libraries (GasSeger PUBLIC pcl_common pcl_kdtree pcl_search pcl_features pcl_io pcl_segmentation pcl_sample_consensus)
# message(STATUS ${PCL_LIBRARIES})
# TODO: Add tests and install targets if needed.
