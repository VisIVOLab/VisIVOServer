cmake_minimum_required(VERSION 3.17)

project(VisIVOServer)

# Set a default build type if none was specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE
    )
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
    )
endif()
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(UTILS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Utils")
set(IMPORTER_SOURCE_DIR "${CMAKE_SOURCE_DIR}/VisIVOImporter")
set(FILTERS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/VisIVOFilters")
set(VSUTILS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/VisIVOUtils")
set(VIEWER_SOURCE_DIR "${CMAKE_SOURCE_DIR}/VisIVOViewer")
set(API_SOURCE_DIR "${CMAKE_SOURCE_DIR}/API")
set(API_LIGHT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/API_LIGHT")

include_directories(${CMAKE_SOURCE_DIR})

option(HAS_CHANGA_IMPORTER "Enable support for Changa importer" OFF) 
option(HAS_HECUBA_IMPORTER "Enable support for Hecuba importer" OFF)
option(INSTALL_API "Enabling VisIVOAPI installation" OFF)
option(API_LIGHT "Install Light version of VisIVOApi" OFF) 

add_subdirectory(${IMPORTER_SOURCE_DIR})
add_subdirectory(${FILTERS_SOURCE_DIR})
add_subdirectory(${VSUTILS_SOURCE_DIR})
add_subdirectory(${VIEWER_SOURCE_DIR})
if(INSTALL_API)
    if(NOT API_LIGHT)
        add_subdirectory(${API_SOURCE_DIR})
    else()
        add_subdirectory(${API_LIGHT_SOURCE_DIR})
    endif()
endif()