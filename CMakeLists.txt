cmake_minimum_required(VERSION 3.18)
project(wacommplusplus)

enable_language(Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)

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

if (USE_MPI AND USE_EMPI)
    message(FATAL_ERROR "You can not use both MPI and EMPI..")
endif()

set(LIBMPI "")
set(LIBEMPI "")
set(LIBPAPI "")
set(LIBGLPK "")
set(LIBICC "")

find_package(MPI)
if(MPI_C_FOUND)
    if (USE_MPI)
        message(STATUS "Using MPI for distributed memory parallelism.")
        add_definitions(-DUSE_MPI)
        include_directories(${MPI_INCLUDE_PATH})
        set(LIBMPI ${MPI_C_LIBRARIES})
    elseif (USE_EMPI)
        find_library(EMPI_LIBRARY empi)
        find_library(PAPI_LIBRARY papi)
        find_library(GLPK_LIBRARY glpk)
        find_library(ICC_LIBRARY icc)

        if (EMPI_LIBRARY AND PAPI_LIBRARY AND GLPK_LIBRARY AND ICC_LIBRARY)
            message(STATUS "Using MPI for distributed memory parallelism with malleability.")
            add_definitions(-DUSE_EMPI)
            include_directories(${MPI_INCLUDE_PATH})
            set(LIBMPI ${MPI_C_LIBRARIES})
        endif()
    endif()
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
	message(STATUS "Use CUDA acceleration.")
	message(STATUS "CUDA Version: ${CUDA_VERSION_STRING}")
	enable_language("CUDA")
	if(${CUDA_VERSION_STRING} VERSION_GREATER_EQUAL "11.1")
      		execute_process(COMMAND ${CMAKE_CUDA_COMPILER} --list-gpu-code RESULT_VARIABLE EXIT_CODE OUTPUT_VARIABLE OUTPUT_VAL)
      		if(EXIT_CODE EQUAL 0)
        		#Remove sm_
        		string(REPLACE "sm_" "" OUTPUT_VAL ${OUTPUT_VAL})
        		#Convert to list
        		string(REPLACE "\n" ";" __CUDA_ARCH_BIN ${OUTPUT_VAL})
        		#Remove last empty entry
        		list(REMOVE_AT __CUDA_ARCH_BIN -1)
      		else()
        		message(FATAL_ERROR "Failed to run NVCC to get list of GPU codes: ${EXIT_CODE}")
      		endif()
    	elseif(${CUDA_VERSION_STRING} VERSION_GREATER_EQUAL "11.0")
      		set(__CUDA_ARCH_BIN "35;37;50;52;53;60;61;62;70;72;75;80")
    	elseif(${CUDA_VERSION_STRING} VERSION_GREATER_EQUAL "10.0")
      		set(__CUDA_ARCH_BIN "30;32;35;37;50;52;53;60;61;62;70;72;75")
    	elseif(${CUDA_VERSION_STRING} VERSION_GREATER_EQUAL "9.1")
      		set(__CUDA_ARCH_BIN "30;32;35;37;50;52;53;60;61;62;70;72")
    	else()
      		set(__CUDA_ARCH_BIN "30;32;35;37;50;52;53;60;61;62;70")
    	endif()
	add_definitions(-DUSE_CUDA)
    
	set(CUDA_ARCH_BIN ${__CUDA_ARCH_BIN} CACHE STRING "Specify 'real' GPU architectures to build binaries for")
    	
    	if(POLICY CMP0104)
		set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_BIN})
    		message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
	else()
		# Generate SASS
    		set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_BIN})
    		# Generate PTX for last architecture
    		list(GET CUDA_ARCH_BIN -1 ver)
    		set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -gencode arch=compute_${ver},code=compute_${ver}")
    		message(STATUS "CMAKE_CUDA_FLAGS: ${CMAKE_CUDA_FLAGS}")
    
    		add_library(pcl_cuda INTERFACE)
    		target_include_directories(pcl_cuda INTERFACE ${CUDA_TOOLKIT_INCLUDE})
	endif()
	#set(CUDA_NVCC_FLAGS "-gencode arch=sm_30")
    	
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
        URL https://zlib.net/zlib-1.3.tar.gz
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
        CONFIGURE_COMMAND ./configure --prefix=${EXTERNAL_INSTALL_LOCATION} CFLAGS=-fPIC CPPFLAGS=-I${EXTERNAL_INSTALL_LOCATION}/include/ LDFLAGS=-L${EXTERNAL_INSTALL_LOCATION}/lib/ --without-libidn --disable-ldap --without-brotli --without-openssl
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
        URL https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/netcdf-c-4.8.1.tar.gz
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
        URL https://downloads.unidata.ucar.edu/netcdf-cxx/4.3.1/netcdf-cxx4-4.3.1.tar.gz
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

add_executable(wacommplusplus main.cpp Particle.cpp Particle.hpp Particles.cpp Particles.hpp Wacomm.cpp Wacomm.hpp Array.h Source.cpp Source.hpp Sources.cpp Sources.hpp OceanModelAdapters/ROMSAdapter.cpp OceanModelAdapters/ROMSAdapter.hpp Config.cpp Config.hpp OceanModelAdapter.cpp OceanModelAdapter.hpp Utils.cpp Utils.hpp WacommPlusPlus.cpp WacommPlusPlus.hpp JulianDate.cpp JulianDate.hpp OceanModelAdapters/WacommAdapter.cpp OceanModelAdapters/WacommAdapter.hpp ${gpu_cu_files} ${gpu_h_files})

# Expicit the dependencies
add_dependencies(zlib szlib)
add_dependencies(curl zlib)
add_dependencies(hdf5 curl)
add_dependencies(netcdf hdf5 )
add_dependencies(netcdfcxx netcdf )
add_dependencies(${PROJECT_NAME} netcdfcxx)
add_dependencies(${PROJECT_NAME} log4cplus)

target_link_libraries(${PROJECT_NAME} ${EMPI_LIBRARY} ${PAPI_LIBRARY} ${GLPK_LIBRARY} ${ICC_LIBRARY} ${LIBOMP} pthread dl nlohmann_json::nlohmann_json ${LIBMPI} ${LIBNETCDFCXX} ${LIBNETCDF} ${LIBCURL} ${LIBHDF5} ${LIBZLIB} ${LIBSZLIB} ${LIBLOG4CPLUS} ${CUDART_LIBRARY} OpenSSL::SSL dl)
