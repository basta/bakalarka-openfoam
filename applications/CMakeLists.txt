cmake_minimum_required(VERSION 3.30)

project(icoHeatSolver)

# Check valid OpenFOAM
if(DEFINED ENV{WM_PROJECT_DIR})
    MESSAGE(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
    message(FATAL_ERROR "The OpenFOAM bashrc is not sourced")
endif(DEFINED ENV{WM_PROJECT_DIR})
set(OpenFOAM_VERSION $ENV{WM_PROJECT_VERSION})
set(OpenFOAM_DIR $ENV{WM_PROJECT_DIR})
set(OpenFOAM_LIB_DIR $ENV{FOAM_LIBBIN})
set(OpenFOAM_SRC $ENV{FOAM_SRC})
# ====================== Some difference between ESI version and Foundation version ==================
# Of course you can change some parameters here, e.g., some pre-settings in ~/.OpenFOAM/prefs.h
set(PATH_LIB_OPENMPI "openmpi-system")  # Foundation version
set(DEFINITIONS_COMPILE "-std=c++11 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedef -Wno-invalid-offsetof -Wno-deprecated-register -Wno-undefined-var-template -O0 -g -DFULLDEBUG -DNoRepository -ftemplate-depth-100 -fPIC")

if(${OpenFOAM_VERSION} MATCHES "v([0-9]*)")       # ESI
    set(PATH_LIB_OPENMPI "sys-openmpi")
    set(DEFINITIONS_COMPILE "-std=c++14 -m64 -pthread -ftrapping-math -DOPENFOAM=2106 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wold-style-cast -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-undefined-var-template -Wno-unknown-warning-option  -O3  -DNoRepository -ftemplate-depth-100 -fPIC -DIMPLEMENT_ACTIVATION -Wl,-execute,-undefined,dynamic_lookup")
endif()
# =====================================================================================================
# Compiling configure
add_definitions("${DEFINITIONS_COMPILE}")
# ======== OS specific setting =============
if(APPLE)
    add_definitions(" -Ddarwin64 ")
endif(APPLE)

# ==========================================

include_directories(.
        ${OpenFOAM_SRC}/OpenFOAM/lnInclude
        ${OpenFOAM_SRC}/OSspecific/POSIX/lnInclude
        ${OpenFOAM_SRC}/finiteVolume/lnInclude
        ${OpenFOAM_SRC}/dynamicFvMesh/lnInclude
        ${OpenFOAM_SRC}/transportModels/incompressible/lnInclude
        ${OpenFOAM_SRC}/TurbulenceModels/incompressible/lnInclude
        ${OpenFOAM_SRC}/thermophysicalModels/radiation/lnInclude
)

link_directories(${OpenFOAM_LIB_DIR} ${OpenFOAM_LIB_DIR}/dummy ${OpenFOAM_LIB_DIR}/${PATH_LIB_OPENMPI})

add_executable(icoHeatSolver
        ./icoHeatExternalForce/buoyantBoussinesqPimpleFoam.C
)