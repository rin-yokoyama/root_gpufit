# CMakeList.txt for the ROOTExamples project
# CMake version 3 is required
cmake_minimum_required (VERSION 3.11 FATAL_ERROR)
project(root_gpufit CXX)

# Add in CMAKE_PREFIX_PATH the installation prefix for ROOT 
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# You need to add COMPONENTS according to what you want to do.
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Minuit)
# Include ROOT cmake macros
include(${ROOT_USE_FILE})

# Set the default install prefix to [source_dir]/install
# (You usually don't want to install this project to /usr/local/)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install CACHE PATH "install dir")

# Set the path to Gpufit library install prefix
# Install from:
# https://github.com/rin-yokoyama/Gpufit.git -b pulse_fit_v1.2.0
set(GPUFIT_PREFIX_PATH ${HOME}/opt/Gpufit CACHE PATH "Path to the Gpufit installation")

function(add_exec cpp_name exec_name)
    add_executable(${exec_name} ${cpp_name})
    target_include_directories(${exec_name} PRIVATE ${GPUFIT_PREFIX_PATH})
    target_link_directories(${exec_name} PRIVATE ${GPUFIT_PREFIX_PATH}/build/Gpufit ${GPUFIT_PREFIX_PATH}/build/Cpufit ${GPUFIT_PREFIX_PATH}/build/Interfaces)
    target_link_libraries(${exec_name} Gpufit Cpufit GpufitInterface ${ROOT_LIBRARIES})
endfunction()

add_exec( pulseFitTestCpuMain.cpp cpu_fit )
add_exec( pulseFitTestGpuMain.cpp gpu_fit )
add_exec( RDataFrameFitMain.cpp rdframe_fit )
