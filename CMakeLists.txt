cmake_minimum_required(VERSION 3.13)
project(wacommplusplus)

enable_language(Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL PGI)
  message(STATUS "Using PGI compilers.")
  add_compile_options(-std=c++17)
endif()

option(DEBUG "Compile logging information in performance-critical classes." OFF)
if(DEBUG)
    add_definitions(-DDEBUG)
endif()

option(USE_MPI "Use MPI for distributed memory parallelism." OFF)
option(USE_OMP "Use OMP for shared memory parallelism." OFF)
option(USE_OPENACC "Use OpenACC acceleration." OFF)
option(USE_CUDA "Use CUDA acceleration." OFF)

set(LIBMPI "")
find_package(MPI)
if(MPI_C_FOUND AND USE_MPI)
    add_definitions(-DUSE_MPI)
    message(STATUS "Using MPI for distributed memory parallelism.")
    include_directories(${MPI_INCLUDE_PATH})
    set(LIBMPI ${MPI_C_LIBRARIES})
endif()

set(LIBOMP "")
find_package(OpenMP)
if(OpenMP_CXX_FOUND AND USE_OMP)
    add_definitions(-DUSE_OMP)
    message(STATUS "Using OMP for shared memory parallelism.")
    set(LIBOMP OpenMP::OpenMP_CXX)
endif()

find_package(OpenACC)
if(OpenACC_FOUND AND USE_OPENACC)
    add_definitions(-DUSE_OPENACC)
    message(STATUS "Use OpenACC acceleration.")
endif()

set(CUDART_LIBRARY "")
set(gpu_h_files "")
set(gpu_cu_files "")
find_package(CUDA)
if(CUDA_FOUND AND USE_CUDA)
    add_definitions(-DUSE_CUDA)
    message(STATUS "Use CUDA acceleration.")
    enable_language("CUDA")
    set(CUDA_NVCC_FLAGS "-arch=sm_30")
    file(GLOB gpu_h_files "${CMAKE_SOURCE_DIR}/cuda/*.h")
    file(GLOB gpu_cu_files "${CMAKE_SOURCE_DIR}/cuda/*.cu")
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()

# Include the support to external projects
include(ExternalProject)

# OpenSSL is needed for building dependences
find_package(OpenSSL REQUIRED)

# Set the external install location
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

set(NETCDFCXXLIBS_PARAMS "")
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(NETCDFCXXLIBS_PARAMS "-lnetcdf -lhdf5_hl -lhdf5 -lm -lz -lcurl")
endif()


# Add to the includes
include_directories(SYSTEM ${EXTERNAL_INSTALL_LOCATION}/include)

# Add to the libraries
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)


include(FetchContent)

FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.7.3)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

ExternalProject_Add(szlib
        URL https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --prefix=${EXTERNAL_INSTALL_LOCATION}
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBSZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libsz.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBSZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libsz.dylib)
endif()

ExternalProject_Add(zlib
        URL https://zlib.net/zlib-1.2.11.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND CFLAGS=-fPIC ./configure --prefix=${EXTERNAL_INSTALL_LOCATION}
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libz.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBZLIB ${EXTERNAL_INSTALL_LOCATION}/lib/libz.dylib)
endif()

ExternalProject_Add(curl
        URL https://curl.haxx.se/download/curl-7.73.0.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_INSTALL_LOCATION} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --without-libidn --disable-ldap
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBCURL ${EXTERNAL_INSTALL_LOCATION}/lib/libcurl.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBCURL ${EXTERNAL_INSTALL_LOCATION}/lib/libcurl.dylib)
endif()

ExternalProject_Add(hdf5
        URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_INSTALL_LOCATION} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --enable-hl --enable-shared --enable-build-mode=production --enable-unsupported --enable-cxx --with-zlib=${EXTERNAL_INSTALL_LOCATION} --with-szlib=${EXTERNAL_INSTALL_LOCATION} --enable-threadsafe --with-pthread
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBHDF5 ${EXTERNAL_INSTALL_LOCATION}/lib/libhdf5_hl.a;${EXTERNAL_INSTALL_LOCATION}/lib/libhdf5.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBHDF5 ${EXTERNAL_INSTALL_LOCATION}/lib/libhdf5.dylib)
endif()

ExternalProject_Add(netcdf
        URL https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --prefix=${EXTERNAL_INSTALL_LOCATION}  --enable-shared --enable-netcdf-4 --enable-dap --enable-byterange --enable-erange-fill
        INSTALL_COMMAND make install
        )
set(LIBNETCDF ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBNETCDF ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf.dylib)
endif()

ExternalProject_Add(netcdfcxx
        URL https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-cxx4-4.3.1.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure LIBS=${NETCDFCXXLIBS_PARAMS} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --prefix=${EXTERNAL_INSTALL_LOCATION}   --enable-shared
        INSTALL_COMMAND make install
        )
set(LIBNETCDFCXX ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf_c++4.a)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBNETCDFCXX ${EXTERNAL_INSTALL_LOCATION}/lib/libnetcdf_c++4.dylib)
endif()


ExternalProject_Add(log4cplus
        URL https://kumisystems.dl.sourceforge.net/project/log4cplus/log4cplus-stable/2.0.5/log4cplus-2.0.5.tar.gz
        TIMEOUT 360
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_INSTALL_LOCATION} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )
set(LIBLOG4CPLUS ${EXTERNAL_INSTALL_LOCATION}/lib/liblog4cplus.so)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(LIBLOG4CPLUS ${EXTERNAL_INSTALL_LOCATION}/lib/liblog4cplus.dylib)
endif()

FetchContent_Declare(cpr GIT_REPOSITORY https://github.com/whoshuu/cpr.git GIT_TAG c8d33915dbd88ad6c92b258869b03aba06587ff9) # the commit hash for 1.5.0
FetchContent_MakeAvailable(cpr)

add_executable(wacommplusplus main.cpp Particle.cpp Particle.hpp Particles.cpp Particles.hpp Wacomm.cpp Wacomm.hpp Array.h Source.cpp Source.hpp Sources.cpp Sources.hpp OceanModelAdapters/ROMSAdapter.cpp OceanModelAdapters/ROMSAdapter.hpp Config.cpp Config.hpp OceanModelAdapter.cpp OceanModelAdapter.hpp Utils.cpp Utils.hpp WacommPlusPlus.cpp WacommPlusPlus.hpp JulianDate.cpp JulianDate.hpp OceanModelAdapters/WacommAdapter.cpp OceanModelAdapters/WacommAdapter.hpp ${gpu_cu_files} ${gpu_h_files})

# Expicit the dependencies
add_dependencies(zlib szlib)
add_dependencies(curl zlib)
add_dependencies(hdf5 curl)
add_dependencies(netcdf hdf5 )
add_dependencies(netcdfcxx netcdf )
add_dependencies(${PROJECT_NAME} netcdfcxx)
add_dependencies(${PROJECT_NAME} log4cplus)

target_link_libraries(${PROJECT_NAME} ${LIBOMP} ${LIBMPI} ${LIBNETCDFCXX} ${LIBNETCDF} ${LIBCURL} ${LIBHDF5} ${LIBZLIB} ${LIBSZLIB} ${LIBLOG4CPLUS} ${CUDART_LIBRARY} pthread dl OpenSSL::SSL nlohmann_json::nlohmann_json cpr::cpr)
