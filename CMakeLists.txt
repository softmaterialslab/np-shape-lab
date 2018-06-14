cmake_minimum_required(VERSION 3.9)
project(Shapes_Simulation)

set(CMAKE_CXX_STANDARD 11)

add_executable(Shapes_Simulation
        control.h
        edge.cpp
        edge.h
        face.cpp
        face.h
        functions.cpp
        functions.h
        interface.cpp
        interface.h
        main.cpp
        newenergies.cpp
        newenergies.h
        newforces.cpp
        newforces.h
        newmd.cpp
        thermostat.h
        utility.h
        vector3d.cpp
        vector3d.h
        vertex.cpp
        vertex.h)

find_package(GSL REQUIRED)
target_link_libraries(Shapes_Simulation GSL::gsl GSL::gslcblas)

find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(Shapes_Simulation ${Boost_LIBRARIES})