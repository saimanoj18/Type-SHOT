cmake_minimum_required(VERSION 2.6 FATAL_ERROR)


project(tro_type_shot)
SET(CMAKE_BUILD_TYPE RelWithDebInfo)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


# Adding PCL library
set(PCL_DIR "/home/sai/workspace/pcl-pcl-1.7.2/build/")
find_package(PCL 1.7.2 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS} include)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


add_executable (shot type_shot/shot.cpp)
target_link_libraries (shot ${PCL_LIBRARIES})


add_executable (type_shot type_shot/type_shot.cpp)
target_link_libraries (type_shot ${PCL_LIBRARIES})


add_executable (type_rops type_rops/type_rops.cpp)
target_link_libraries (type_rops ${PCL_LIBRARIES})


add_executable (type_fpfh type_fpfh/type_fpfh.cpp)
target_link_libraries (type_fpfh ${PCL_LIBRARIES})





