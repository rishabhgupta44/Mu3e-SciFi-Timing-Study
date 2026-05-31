cmake_minimum_required(VERSION 3.16)
project(Mu3eSciFiTimingStudy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED ui_all vis_all)

list(GET Geant4_INCLUDE_DIRS 0 Geant4_INCLUDE_DIR)
get_filename_component(_geant4_include_dir "${Geant4_INCLUDE_DIR}" DIRECTORY)
get_filename_component(_geant4_prefix_dir "${_geant4_include_dir}" DIRECTORY)
set(_geant4_lib_dir "${_geant4_prefix_dir}/lib")

set(CMAKE_BUILD_RPATH "${_geant4_lib_dir}")
set(CMAKE_INSTALL_RPATH "${_geant4_lib_dir}")

add_executable(scifi_sim src/main.cc)

target_include_directories(scifi_sim PRIVATE ${Geant4_INCLUDE_DIRS})
target_link_libraries(scifi_sim ${Geant4_LIBRARIES})
