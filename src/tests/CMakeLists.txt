
# ***********************************************************************
# Copyright (c) 2015 Advanced Micro Devices, Inc. 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions 
# are met:
# 
# 1. Redistributions of source code must retain the above copyright 
# notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# ***********************************************************************

set( clProbDist.CTest.Source
  ctest/cli.c
  ctest/util.c
  ctest/exponential_dispatch.c
  ctest/exponential_checks.c
  ctest/normal_dispatch.c
  ctest/normal_checks.c
  ctest/lognormal_dispatch.c
  ctest/lognormal_checks.c
  ctest/gamma_dispatch.c
  ctest/gamma_checks.c
  ctest/poisson_dispatch.c
  ctest/poisson_checks.c
  ${PROJECT_SOURCE_DIR}/client/common.c
)
set( clProbDist.CTest.Headers
  ${PROJECT_SOURCE_DIR}/include/clProbDist/clProbDist.h
  ${PROJECT_SOURCE_DIR}/include/clProbDist/exponential.h
  ${PROJECT_SOURCE_DIR}/include/clProbDist/gamma.h
  ${PROJECT_SOURCE_DIR}/include/clProbDist/normal.h
  ${PROJECT_SOURCE_DIR}/include/clProbDist/lognormal.h
  ${PROJECT_SOURCE_DIR}/include/clProbDist/poisson.h
  ${PROJECT_SOURCE_DIR}/client/common.h
  ctest/util.h
  ctest/mangle.h
  ctest/checks.h
  ctest/dispatch.c.h
  ctest/checks.c.h
)

set( clProbDist.CTest.Files ${clProbDist.CTest.Source} ${clProbDist.CTest.Headers} )

include_directories(
  ${OPENCL_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
)

set( Client.Source ${clProbDist.CTest.Source} )

if( MSVC )
	if( MSVC_VERSION LESS 1800 )
		# Use C++ with Microsoft compiler
		SET_SOURCE_FILES_PROPERTIES( ${Client.Source} PROPERTIES LANGUAGE CXX)
	endif ()
endif( )
		
set( DL_LIB "" )
if( WIN32 )
        add_definitions( "/D_CONSOLE" )
elseif( APPLE )
        set( CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}" )   
else( )
        # To use the dlopen() and dlclose() functions, we should link with libdl
        set( DL_LIB "-ldl" )
endif( )

if( CMAKE_COMPILER_IS_GNUCC )
    set( MATH_LIB "-lm" )
endif()

if( NOT BUILD_RUNTIME )
  message( ERROR "Building tests requires building runtime." )
endif( )

add_executable( CTest   ${clProbDist.CTest.Files} )

include_directories(   CTest               ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
target_link_libraries( CTest    clProbDist ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

set_target_properties( CTest    PROPERTIES VERSION ${CLPROBDIST_VERSION} )
set_target_properties( CTest    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )


# CPack configuration; include the executable into the package
install( TARGETS CTest
        RUNTIME DESTINATION bin${SUFFIX_BIN}
        )

