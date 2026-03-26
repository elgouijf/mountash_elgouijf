cmake_minimum_required(VERSION 3.10)
project(SimuParticules)

set(CMAKE_CXX_STANDARD 17)

add_executable(Particule Particule.cxx)

add_custom_target(run
    COMMAND ./Particule > ${CMAKE_SOURCE_DIR}/systeme_solaire.txt
    COMMAND ${CMAKE_SOURCE_DIR}/script_visual.py
    DEPENDS Particule
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_subdirectory(doc)