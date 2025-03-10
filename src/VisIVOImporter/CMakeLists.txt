cmake_minimum_required(VERSION 3.15)

project(VisIVOImporter CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(CFITSIO REQUIRED)
find_package(HDF5 REQUIRED)
find_package(CURL REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

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
    gl2ps
    glew
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

find_library(TIRPC_LIBRARIES
    NAMES tirpc
)

set(IMPORTER_SOURCE
    Importers/abstractsource.cpp
    Importers/abstractsource.h
    Importers/asciisource.cpp
    Importers/asciisource.h
    Importers/binsource.cpp
    Importers/binsource.h
    Importers/changasource.cpp
    Importers/changasource.h
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
)

find_package(HDF5 REQUIRED)

add_executable(VisIVOImporter ${IMPORTER_SOURCE} ${UTILS_SOURCE})

target_include_directories(VisIVOImporter PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/Importers
    ${UTILS_SOURCE_DIR}
)

target_link_libraries(VisIVOImporter PRIVATE
    CFITSIO::CFITSIO
    CURL::libcurl
    ${MPI_C_LIBRARIES}
    OpenMP::OpenMP_CXX
    ${VTK_LIBRARIES}
    stdc++fs
    ${TIRPC_LIBRARIES}
)

install(TARGETS VisIVOImporter DESTINATION bin)

target_compile_definitions(VisIVOImporter PRIVATE LIGHT)

set(CXXFLAGS  -fopenmp)
set(CMAKE_CXX_FLAGS "${CXXFLAGS}")
target_link_libraries(VisIVOImporter PUBLIC MPI::MPI_CXX)
