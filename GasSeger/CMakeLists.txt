# CMakeList.txt : CMake project for GasSeger, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (GasSeger "GasSeger.cpp" "Calc.cpp" "GasSeger.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET GasSeger PROPERTY CXX_STANDARD 20)
endif()

find_package(PCL REQUIRED )

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
link_libraries(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

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
