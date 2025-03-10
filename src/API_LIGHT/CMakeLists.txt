cmake_minimum_required(VERSION 3.15)

project(VisIVOApi
    VERSION 2.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_INCLUDE_CURRENT_DIR True)


find_package(Boost REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(HDF5 REQUIRED)
find_package(CURL REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_library(TIRPC_LIBRARIES
    NAMES tirpc
)


set(API_HEADERS
    visivo.h
    visivodef.h
    visivoserver.h
)

set(API_SOURCE
    visivoimp.cpp
    visivofil.cpp
    visivosetLight.cpp
)

set(IMPORTER_SOURCE
    ${IMPORTER_SOURCE_DIR}/Importers/abstractsource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/abstractsource.h
    ${IMPORTER_SOURCE_DIR}/Importers/asciisource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/asciisource.h
    ${IMPORTER_SOURCE_DIR}/Importers/binsource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/binsource.h
    ${IMPORTER_SOURCE_DIR}/Importers/changasource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/changasource.h
    ${IMPORTER_SOURCE_DIR}/Importers/csvsource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/csvsource.h
    ${IMPORTER_SOURCE_DIR}/Importers/fitsimagesource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/fitsimagesource.h
    ${IMPORTER_SOURCE_DIR}/Importers/fitstablesource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/fitstablesource.h
    ${IMPORTER_SOURCE_DIR}/Importers/flysource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/flysource.h
    ${IMPORTER_SOURCE_DIR}/Importers/gadgetsource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/gadgetsource.h
    ${IMPORTER_SOURCE_DIR}/Importers/hdf5source.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/hdf5source.h
    ${IMPORTER_SOURCE_DIR}/Importers/muportalsource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/muportalsource.h
    ${IMPORTER_SOURCE_DIR}/Importers/ramsessource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/ramsessource.h
    ${IMPORTER_SOURCE_DIR}/Importers/rawgridsource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/rawgridsource.h
    ${IMPORTER_SOURCE_DIR}/Importers/rawpointssource.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/rawpointssource.h
    ${IMPORTER_SOURCE_DIR}/Importers/vosourcenew.cpp
    ${IMPORTER_SOURCE_DIR}/Importers/vosourcenew.h
)

set(FILTERS_SOURCE
    ${FILTERS_SOURCE_DIR}/Filters/fpoptimizer.cpp
    ${FILTERS_SOURCE_DIR}/Filters/startFilter.cpp
    ${FILTERS_SOURCE_DIR}/Filters/startFilter.h
    ${FILTERS_SOURCE_DIR}/Filters/vsaddidop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsaddidop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsappend.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsappend.h
    ${FILTERS_SOURCE_DIR}/Filters/vschangecolnameop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vschangecolnameop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsclustercellop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsclustercellop.h
    ${FILTERS_SOURCE_DIR}/Filters/vscoarsevolumeop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vscoarsevolumeop.h
    ${FILTERS_SOURCE_DIR}/Filters/vscutop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vscutop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsdecimatorop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsdecimatorop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsexampleop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsexampleop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsextractlistrowsop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsextractlistrowsop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsextractsubvolumeop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsextractsubvolumeop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsgrid2pointdistr.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsgrid2pointdistr.h
    ${FILTERS_SOURCE_DIR}/Filters/vsincludeop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsincludeop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsinterpolateop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsinterpolateop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsLine.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsLine.h
    ${FILTERS_SOURCE_DIR}/Filters/vsmathop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsmathop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsmergeop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsmergeop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsmoduleop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsmoduleop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsmrcampos.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsmrcampos.h
    ${FILTERS_SOURCE_DIR}/Filters/vsmuportalop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsmuportalop.h
    ${FILTERS_SOURCE_DIR}/Filters/vspointdistributeop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vspointdistributeop.h
    ${FILTERS_SOURCE_DIR}/Filters/vspointpropertyop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vspointpropertyop.h
    ${FILTERS_SOURCE_DIR}/Filters/vspolarop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vspolarop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsprinttableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsprinttableop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsrandomizertableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsrandomizertableop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsselectcolumnsop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsselectcolumnsop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsselectfieldop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsselectfieldop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsselectvol.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsselectvol.h
    ${FILTERS_SOURCE_DIR}/Filters/vsshowtableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsshowtableop.h
    ${FILTERS_SOURCE_DIR}/Filters/vssigmacontoursop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vssigmacontoursop.h
    ${FILTERS_SOURCE_DIR}/Filters/vssphereop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vssphereop.h
    ${FILTERS_SOURCE_DIR}/Filters/vssplittableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vssplittableop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsstatisticop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsstatisticop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsswapop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsswapop.h
    ${FILTERS_SOURCE_DIR}/Filters/vstableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vstableop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsVector.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsVector.h
    ${FILTERS_SOURCE_DIR}/Filters/vsvisualop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsvisualop.h
    ${FILTERS_SOURCE_DIR}/Filters/vsvollimit.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsvollimit.h
    ${FILTERS_SOURCE_DIR}/Filters/vsVoxel.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vsVoxel.h
    ${FILTERS_SOURCE_DIR}/Filters/vswrvotableop.cpp
    ${FILTERS_SOURCE_DIR}/Filters/vswrvotableop.h
)

set(VSUTILS_SOURCE
    ${VSUTILS_SOURCE_DIR}/vscreategenericslices.cpp
    ${VSUTILS_SOURCE_DIR}/vscreategenericslices.h
    ${VSUTILS_SOURCE_DIR}/vscreatepath.cpp
    ${VSUTILS_SOURCE_DIR}/vscreatepath.h
    ${VSUTILS_SOURCE_DIR}/vscreateslices.cpp
    ${VSUTILS_SOURCE_DIR}/vscreateslices.h
    ${VSUTILS_SOURCE_DIR}/vsloadhistory.cpp
    ${VSUTILS_SOURCE_DIR}/vsloadhistory.h
    ${VSUTILS_SOURCE_DIR}/vstextcol.cpp
    ${VSUTILS_SOURCE_DIR}/vstextcol.h
    ${VSUTILS_SOURCE_DIR}/vsutils.cpp
    ${VSUTILS_SOURCE_DIR}/vsutils.h
    ${VSUTILS_SOURCE_DIR}/vuparametersparser.cpp
    ${VSUTILS_SOURCE_DIR}/vuparametersparser.h
)

set(UTILS_SOURCE
    ${UTILS_SOURCE_DIR}/commandline.cpp
    ${UTILS_SOURCE_DIR}/commandline.h
    ${UTILS_SOURCE_DIR}/fparser.cpp
    ${UTILS_SOURCE_DIR}/fparser.h
    ${UTILS_SOURCE_DIR}/historyXmlWriter.cc
    ${UTILS_SOURCE_DIR}/historyXmlWriter.h
    ${UTILS_SOURCE_DIR}/parametersparser.cpp
    ${UTILS_SOURCE_DIR}/parametersparser.h
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
)

add_library(VisIVOApi SHARED
    ${API_SOURCE}
    ${API_HEADERS}
    ${IMPORTER_SOURCE}
    ${FILTERS_SOURCE}
    ${VSUTILS_SOURCE}
    ${UTILS_SOURCE}
)

set_target_properties(VisIVOApi PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_definitions(VisIVOApi PRIVATE MAC LIGHT)

target_include_directories(VisIVOApi PRIVATE
    ${IMPORTER_SOURCE_DIR}/Importers
    ${FILTERS_SOURCE_DIR}/Filters
    ${VSUTILS_SOURCE_DIR}
    ${UTILS_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(VisIVOApi PRIVATE
    CFITSIO::CFITSIO
    CURL::libcurl
    ${MPI_C_LIBRARIES}
    OpenMP::OpenMP_CXX
    ${VTK_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${TIRPC_LIBRARIES}
)

install(TARGETS VisIVOApi DESTINATION lib)
install(FILES ${API_HEADERS} DESTINATION include)

set(CXXFLAGS  -fopenmp)
set(CMAKE_CXX_FLAGS "${CXXFLAGS}")
target_link_libraries(VisIVOApi PUBLIC MPI::MPI_CXX)
