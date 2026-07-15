cmake_minimum_required(VERSION 3.15)

project(VisIVOFilter CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Boost REQUIRED)

set(FILTERS_SOURCE
    Filters/fpoptimizer.cpp
    Filters/mainFilter.cpp
    Filters/startFilter.cpp
    Filters/startFilter.h
    Filters/vsaddidop.cpp
    Filters/vsaddidop.h
    Filters/vsaddrowsop.cpp
    Filters/vsaddrowsop.h
    Filters/vsappend.cpp
    Filters/vsappend.h
    Filters/vschangecolnameop.cpp
    Filters/vschangecolnameop.h
    Filters/vsclustercellop.cpp
    Filters/vsclustercellop.h
    Filters/vscoarsevolumeop.cpp
    Filters/vscoarsevolumeop.h
    Filters/vscutop.cpp
    Filters/vscutop.h
    Filters/vsdecimatorop.cpp
    Filters/vsdecimatorop.h
    Filters/vsexampleop.cpp
    Filters/vsexampleop.h
    Filters/vsextractlistrowsop.cpp
    Filters/vsextractlistrowsop.h
    Filters/vsextractsubvolumeop.cpp
    Filters/vsextractsubvolumeop.h
    Filters/vsgrid2pointdistr.cpp
    Filters/vsgrid2pointdistr.h
    Filters/vsincludeop.cpp
    Filters/vsincludeop.h
    Filters/vsinterpolateop.cpp
    Filters/vsinterpolateop.h
    Filters/vsLine.cpp
    Filters/vsLine.h
    Filters/vsmathop.cpp
    Filters/vsmathop.h
    Filters/vsmergeop.cpp
    Filters/vsmergeop.h
    Filters/vsmoduleop.cpp
    Filters/vsmoduleop.h
    Filters/vsmrcampos.cpp
    Filters/vsmrcampos.h
    Filters/vsmuportalop.cpp
    Filters/vsmuportalop.h
    Filters/vspointdistributeop.cpp
    Filters/vspointdistributeop.h
    Filters/vspointpropertyop.cpp
    Filters/vspointpropertyop.h
    Filters/vspolarop.cpp
    Filters/vspolarop.h
    Filters/vsprinttableop.cpp
    Filters/vsprinttableop.h
    Filters/vsrandomizertableop.cpp
    Filters/vsrandomizertableop.h
    Filters/vsselectcolumnsop.cpp
    Filters/vsselectcolumnsop.h
    Filters/vsselectfieldop.cpp
    Filters/vsselectfieldop.h
    Filters/vsselectvol.cpp
    Filters/vsselectvol.h
    Filters/vsshowtableop.cpp
    Filters/vsshowtableop.h
    Filters/vssigmacontoursop.cpp
    Filters/vssigmacontoursop.h
    Filters/vssphereop.cpp
    Filters/vssphereop.h
    Filters/vssplittableop.cpp
    Filters/vssplittableop.h
    Filters/vsstatisticop.cpp
    Filters/vsstatisticop.h
    Filters/vsswapop.cpp
    Filters/vsswapop.h
    Filters/vstableop.cpp
    Filters/vstableop.h
    Filters/vsVector.cpp
    Filters/vsVector.h
    Filters/vsvisualop.cpp
    Filters/vsvisualop.h
    Filters/vsvollimit.cpp
    Filters/vsvollimit.h
    Filters/vsVoxel.cpp
    Filters/vsVoxel.h
    Filters/vswrvotableop.cpp
    Filters/vswrvotableop.h
)

set(UTILS_SOURCE
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
    ${UTILS_SOURCE_DIR}/util.cpp
    ${UTILS_SOURCE_DIR}/util.h
    ${UTILS_SOURCE_DIR}/vsobject.cpp
    ${UTILS_SOURCE_DIR}/vsobject.h
    ${UTILS_SOURCE_DIR}/vstable.cpp
    ${UTILS_SOURCE_DIR}/vstable.h
)

add_executable(VisIVOFilter ${FILTERS_SOURCE} ${UTILS_SOURCE})

install(TARGETS VisIVOFilter DESTINATION bin)

target_include_directories(VisIVOFilter PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/Filters
    ${UTILS_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
)
