cmake_minimum_required(VERSION 3.15)

project(VisIVOViewer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(CURL REQUIRED)
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

message(STATUS "VisIVOViewer: VTK version " ${VTK_VERSION})

set(VIEWER_SOURCE
    Viewer/isosurfacepipe.cpp
    Viewer/isosurfacepipe.h
    Viewer/mainViewer.cpp
    Viewer/pipe.cpp
    Viewer/pipe.h
    Viewer/pointspipe.cpp
    Viewer/pointspipe.h
    Viewer/pointssmoothpipe.cpp
    Viewer/pointssmoothpipe.h
    Viewer/slicerpipe.cpp
    Viewer/slicerpipe.h
    Viewer/vectorpipe.cpp
    Viewer/vectorpipe.h
    Viewer/volumepipe.cpp
    Viewer/volumepipe.h
    Viewer/vtkimagepipe.cpp
    Viewer/vtkimagepipe.h
)

set(FILTERS_SOURCE
    ${FILTERS_SOURCE_DIR}/Filters/vsstatisticop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsstatisticop.h
    ${FILTERS_SOURCE_DIR}/Filters/vstableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vstableop.h
)

set(UTILS_SOURCE
    ${UTILS_SOURCE_DIR}/color.cpp
    ${UTILS_SOURCE_DIR}/color.h
    ${UTILS_SOURCE_DIR}/extendedglyph3d.cpp
    ${UTILS_SOURCE_DIR}/extendedglyph3d.h
    ${UTILS_SOURCE_DIR}/historyXmlWriter.cc
    ${UTILS_SOURCE_DIR}/historyXmlWriter.h
    ${UTILS_SOURCE_DIR}/luteditor.cpp
    ${UTILS_SOURCE_DIR}/luteditor.h
    ${UTILS_SOURCE_DIR}/optionssetter.cpp
    ${UTILS_SOURCE_DIR}/optionssetter.h
    ${UTILS_SOURCE_DIR}/tinystr.cpp
    ${UTILS_SOURCE_DIR}/tinystr.h
    ${UTILS_SOURCE_DIR}/tinyxml.cpp
    ${UTILS_SOURCE_DIR}/tinyxml.h
    ${UTILS_SOURCE_DIR}/tinyxmlerror.cpp
    ${UTILS_SOURCE_DIR}/tinyxmlparser.cpp
    ${UTILS_SOURCE_DIR}/visivoutils.cpp
    ${UTILS_SOURCE_DIR}/visivoutils.h
    ${UTILS_SOURCE_DIR}/vsobject.cpp
    ${UTILS_SOURCE_DIR}/vsobject.h
    ${UTILS_SOURCE_DIR}/vstable.cpp
    ${UTILS_SOURCE_DIR}/vstable.h
    ${UTILS_SOURCE_DIR}/vstablemem.cpp
    ${UTILS_SOURCE_DIR}/vstablemem.h
)

add_executable(VisIVOViewer ${VIEWER_SOURCE} ${FILTERS_SOURCE} ${UTILS_SOURCE})

install(TARGETS VisIVOViewer DESTINATION bin)

target_compile_definitions(VisIVOViewer PRIVATE LIGHT)

target_include_directories(VisIVOViewer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Viewer
    ${FILTERS_SOURCE_DIR}/Filters
    ${UTILS_SOURCE_DIR}
    ${VSUTILS_SOURCE_DIR}
)

target_link_libraries(VisIVOViewer PRIVATE
    ${VTK_LIBRARIES}
    CURL::libcurl
)

vtk_module_autoinit(
    TARGETS VisIVOViewer
    MODULES ${VTK_LIBRARIES}
)
