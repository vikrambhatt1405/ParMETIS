cmake_minimum_required(VERSION 3.15)
project(ParMetis)
set(CMAKE_CXX_STANDARD 14)
message(${CMAKE_VERSION})
include_directories(.)
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
add_executable(ParMetis
        adjtocsr.cpp
        CMakeLists.txt
        main.cpp)
target_link_libraries(ParMetis ${MPI_LIBRARIES} libparmetis.a libmetis.a)
if(MPI_COMPILE_FLAGS)
    set_target_properties(ParMetis PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(ParMetis PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()