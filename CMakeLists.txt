cmake_minimum_required(VERSION 3.12)
project(Sofa.Fenics LANGUAGES CXX)

set(SOFA_FENICS_SRC "src/sofa/fenics")

find_package(Sofa.Core REQUIRED)

set(HEADER_FILES
    ${SOFA_FENICS_SRC}/config.h.in
    ${SOFA_FENICS_SRC}/forcefield/HyperElasticForceField.h
)

set(SOURCE_FILES
    ${SOFA_FENICS_SRC}/init.cpp
    ${SOFA_FENICS_SRC}/forcefield/HyperElasticForceField.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${Sofa_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "${PROJECT_NAME}"
    RELOCATABLE "plugins"
)
