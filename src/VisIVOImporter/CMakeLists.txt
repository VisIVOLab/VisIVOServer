cmake_minimum_required(VERSION 3.15)

project(VisIVOImporter CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(CFITSIO REQUIRED)
find_package(HDF5 REQUIRED)
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(MPI REQUIRED)

if(APPLE)
    message(STATUS "Configuring OpenMP for macOS (Apple Clang).")

    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    set(OpenMP_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
    set(OpenMP_FOUND TRUE)
else()
    find_package(OpenMP REQUIRED)
endif()

find_package(VTK 9.1.0 REQUIRED COMPONENTS
    ChartsCore
    CommonColor
    CommonComputationalGeometry
    CommonCore
    CommonDataModel
    CommonExecutionModel
    CommonMath
    CommonMisc
    CommonSystem
    CommonTransforms
    DICOMParser
    DomainsChemistry
    doubleconversion
    exodusII
    expat
    FiltersAMR
    FiltersCore
    FiltersExtraction
    FiltersFlowPaths
    FiltersGeneral
    FiltersGeneric
    FiltersGeometry
    FiltersHybrid
    FiltersHyperTree
    FiltersImaging
    FiltersModeling
    FiltersParallel
    FiltersParallelImaging
    FiltersProgrammable
    FiltersSelection
    FiltersSMP
    FiltersSources
    FiltersStatistics
    FiltersTexture
    FiltersVerdict
    freetype
    GeovisCore
    glew
    gl2ps
    ImagingColor
    ImagingCore
    ImagingFourier
    ImagingGeneral
    ImagingHybrid
    ImagingMath
    ImagingMorphological
    ImagingSources
    ImagingStatistics
    ImagingStencil
    InfovisCore
    InfovisLayout
    InteractionImage
    InteractionStyle
    InteractionWidgets
    IOAMR
    IOCore
    IOEnSight
    IOExodus
    IOExport
    IOGeometry
    IOImage
    IOImport
    IOInfovis
    IOLegacy
    IOLSDyna
    IOMINC
    IOMovie
    IOParallel
    IOParallelXML
    IOPLY
    IOSQL
    IOVideo
    IOXML
    IOXMLParser
    jsoncpp
    loguru
    lz4
    lzma
    metaio
    ogg
    ParallelCore
    ParallelDIY
    RenderingAnnotation
    RenderingContext2D
    RenderingCore
    RenderingFreeType
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingImage
    RenderingLabel
    RenderingLOD
    RenderingOpenGL2
    RenderingUI
    RenderingVolume
    RenderingVolumeOpenGL2
    sqlite
    vtksys
    verdict
    ViewsContext2D
    ViewsCore
    ViewsInfovis
    zlib
)

set(IMPORTER_SOURCE
    Importers/abstractsource.cpp
    Importers/abstractsource.h
    Importers/asciisource.cpp
    Importers/asciisource.h
    Importers/binsource.cpp
    Importers/binsource.h
    Importers/csvsource.cpp
    Importers/csvsource.h
    Importers/fitsimagesource.cpp
    Importers/fitsimagesource.h
    Importers/fitstablesource.cpp
    Importers/fitstablesource.h
    Importers/flysource.cpp
    Importers/flysource.h
    Importers/gadgetsource.cpp
    Importers/gadgetsource.h
    Importers/hdf5source.cpp
    Importers/hdf5source.h
    Importers/mainImporter.cpp
    Importers/muportalsource.cpp
    Importers/muportalsource.h
    Importers/ramsessource.cpp
    Importers/ramsessource.h
    Importers/rawgridsource.cpp
    Importers/rawgridsource.h
    Importers/rawpointssource.cpp
    Importers/rawpointssource.h
    Importers/vosourcenew.cpp
    Importers/vosourcenew.h
    Importers/vtksource.cpp
    Importers/vtksource.h
)

if(HAS_CHANGA_IMPORTER)
    list(APPEND IMPORTER_SOURCE
        Importers/changasource.cpp
        Importers/changasource.h
    )
endif()

if(HAS_HECUBA_IMPORTER)
    list(APPEND IMPORTER_SOURCE
        Importers/hecubasource.cpp
        Importers/hecubasource.h
    )
endif()

set(UTILS_SOURCE
    ${UTILS_SOURCE_DIR}/commandline.cpp
    ${UTILS_SOURCE_DIR}/commandline.h
    ${UTILS_SOURCE_DIR}/historyXmlWriter.cc
    ${UTILS_SOURCE_DIR}/historyXmlWriter.h
    ${UTILS_SOURCE_DIR}/tinystr.cpp
    ${UTILS_SOURCE_DIR}/tinystr.h
    ${UTILS_SOURCE_DIR}/tinyxml.cpp
    ${UTILS_SOURCE_DIR}/tinyxml.h
    ${UTILS_SOURCE_DIR}/tinyxmlerror.cpp
    ${UTILS_SOURCE_DIR}/tinyxmlparser.cpp
    ${UTILS_SOURCE_DIR}/visivoutils.cpp
    ${UTILS_SOURCE_DIR}/visivoutils.h
    ${UTILS_SOURCE_DIR}/vstable.h
    ${UTILS_SOURCE_DIR}/vstable.cpp
    ${UTILS_SOURCE_DIR}/vstablemem.h
    ${UTILS_SOURCE_DIR}/vstablemem.cpp
    ${UTILS_SOURCE_DIR}/vsobject.h
    ${UTILS_SOURCE_DIR}/vsobject.cpp
)

add_executable(VisIVOImporter ${IMPORTER_SOURCE} ${UTILS_SOURCE})

target_include_directories(VisIVOImporter PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/Importers
    ${UTILS_SOURCE_DIR}
)

if(APPLE AND OpenMP_INCLUDE_DIR)
    target_include_directories(VisIVOImporter PRIVATE ${OpenMP_INCLUDE_DIR})
endif()

target_link_libraries(VisIVOImporter PRIVATE
    CFITSIO::CFITSIO
    CURL::libcurl
    ZLIB::ZLIB
    ${MPI_C_LIBRARIES}
    ${VTK_LIBRARIES}
)

if(APPLE)
    target_link_libraries(VisIVOImporter PRIVATE ${OpenMP_omp_LIBRARY})
    target_compile_options(VisIVOImporter PRIVATE "-Xpreprocessor" "-fopenmp")
else()
    target_link_libraries(VisIVOImporter PRIVATE OpenMP::OpenMP_CXX)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT APPLE)
    target_link_libraries(VisIVOImporter PRIVATE stdc++fs)
endif()

if(HAS_CHANGA_IMPORTER)
    find_library(TIRPC_LIBRARIES
        NAMES tirpc
    )
    target_link_libraries(VisIVOImporter PRIVATE ${TIRPC_LIBRARIES})
    target_compile_definitions(VisIVOImporter PRIVATE HAS_CHANGA_IMPORTER=1)
endif()

if(HAS_HECUBA_IMPORTER)
    set(HECUBA_ROOT "/usr/local" CACHE PATH "Path to Hecuba installation")

    find_library(HECUBA_HFETCH_LIBRARY hfetch PATHS ${HECUBA_ROOT}/lib)
    if(NOT HECUBA_HFETCH_LIBRARY)
        message(FATAL_ERROR "Hecuba hfetch library not found in ${HECUBA_ROOT}/lib")
    endif()
    target_include_directories(VisIVOImporter PRIVATE ${HECUBA_ROOT}/include ${HECUBA_ROOT}/include/hecuba)
    target_link_libraries(VisIVOImporter PRIVATE ${HECUBA_HFETCH_LIBRARY}) 
    target_compile_definitions(VisIVOImporter PRIVATE HAS_HECUBA_IMPORTER=1)
endif()

install(TARGETS VisIVOImporter DESTINATION bin)

target_compile_definitions(VisIVOImporter PRIVATE LIGHT)

target_link_libraries(VisIVOImporter PUBLIC MPI::MPI_CXX)

if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(VisIVOImporter PRIVATE "-stdlib=libc++")
    target_link_options(VisIVOImporter PRIVATE "-stdlib=libc++")
endif()
