cmake_minimum_required(VERSION 3.5)

project(fea LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release)

find_package(VTK REQUIRED)
set(VTKFIG_INCLUDE_DIRS "@CONF_INCLUDE_DIRS@;@VTK_INCLUDE_DIRS@;@VTK_PREFIX_PATH@/include/vtk-@VTK_MAJOR_VERSION@.@VTK_MINOR_VERSION@")

find_package(PCL 1.12 REQUIRED)

INCLUDE_DIRECTORIES(
    ${PCL_INCLUDE_DIRS}
    ${VTKFIG_INCLUDE_DIRS}
)

link_directories(
    ${PCL_LIBRARY_DIRS}
    ${PCL_LIBRARY_DIRS}/..
)

add_definitions(
    ${PCL_DEFINITIONS}
)


add_library(fea ${G2O_LIB_TYPE}
    pos.cpp
    pos.hpp
    fem.cpp
    fem.hpp
    fea.cpp
    fea.hpp
)

target_include_directories(fea PUBLIC
  "$<BUILD_INTERFACE:${g2o_SOURCE_DIR};${PROJECT_BINARY_DIR}>"
  $<INSTALL_INTERFACE:include/g2o/fea>
)

set_target_properties(fea PROPERTIES OUTPUT_NAME ${LIB_PREFIX}fea)
set_target_properties(fea PROPERTIES
  VERSION ${G2O_LIB_VERSION}
  SOVERSION ${G2O_LIB_SOVERSION})
target_link_libraries(fea PUBLIC 
    ${G2O_EIGEN3_EIGEN_TARGET}
    ${VTK_LIBRARIES}
    ${PCL_LIBRARIES})
target_compile_features(fea PUBLIC cxx_std_17)

install(TARGETS fea
  EXPORT ${G2O_TARGETS_EXPORT_NAME}
  RUNTIME DESTINATION ${RUNTIME_DESTINATION}
  LIBRARY DESTINATION ${LIBRARY_DESTINATION}
  ARCHIVE DESTINATION ${ARCHIVE_DESTINATION}
  INCLUDES DESTINATION ${INCLUDES_DESTINATION}
)

file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

install(FILES ${headers} DESTINATION ${INCLUDES_INSTALL_DIR}/fea)