# ########################################################################



set( Common.Headers common.h
    ../include/clProbDist/clProbDist.h
    ../include/clProbDist/poisson.h
    ../include/clProbDist/gamma.h
    ../include/clProbDist/normal.h
    ../include/clProbDist/exponential.h
    )
set( Common.Source common.c )

# Self-Contained Example
set( SelfContained.Source    selfcontained.c )
set( SelfContained.Files     ${SelfContained.Source} )

# Docs Tutorial
set( DocsTutorial1.Source    DocsTutorial/example1.c ${Common.Source} )
set( DocsTutorial1.Headers   ${Common.Headers} )
set( DocsTutorial1.Files     ${DocsTutorial1.Source} ${DocsTutorial1.Headers} )
set( DocsTutorial2.Source    DocsTutorial/example2.c ${Common.Source} )
set( DocsTutorial2.Headers   ${Common.Headers} )
set( DocsTutorial2.Files     ${DocsTutorial2.Source} ${DocsTutorial2.Headers} )
set( DocsTutorial3.Source    DocsTutorial/example3.c ${Common.Source} )
set( DocsTutorial3.Headers   ${Common.Headers} )
set( DocsTutorial3.Files     ${DocsTutorial3.Source} ${DocsTutorial3.Headers} )
set( DocsTutorial4.Source    DocsTutorial/example4.c ${Common.Source} )
set( DocsTutorial4.Headers   ${Common.Headers} )
set( DocsTutorial4.Files     ${DocsTutorial4.Source} ${DocsTutorial4.Headers}
                             DocsTutorial/example4_kernel.cl )

# Inventory
set( Inventory.Source
        Inventory/inventory.c
        Inventory/Policies.c
        Inventory/SimulateRuns.c
        ${Common.Source} )
set( Inventory.Headers
        Inventory/Types.h
        Inventory/Policies.h
        Inventory/SimulateRuns.h
        ${Common.Headers} )
set( Inventory.Files ${Inventory.Source} ${Inventory.Headers}
        Inventory/InventoryKernels.cl )

# Asian Option
set( AsianOption.Source
        AsianOption/AsianOption.c
        AsianOption/Simulate_Runs.c
        ${Common.Source} )
set( AsianOption.Headers
        AsianOption/Types.h
        AsianOption/Simulate_Runs.h
        ${Common.Headers} )
set( AsianOption.Files ${AsianOption.Source} ${AsianOption.Headers}
        AsianOption/AsianOptionKernel.cl )

# Performance
set( Performance.Source
    Performance/performance.c
    Performance/SimulateRun.c
    Performance/Options.c
        ${Common.Source} )
    set( Performance.Headers
        Performance/Options.h
        Performance/SimulateRun.h
        ${Common.Headers} )
    set( Performance.Files ${Performance.Source} ${Performance.Headers}
        Performance/PerformanceKernels.cl )

set( Client.Source ${Inventory.Source} )

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

# Include standard OpenCL headers
include_directories( ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )

# CPack configuration; include the executable into the package

add_executable( SelfContained   ${SelfContained.Files} )

include_directories(   SelfContained               ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
target_link_libraries( SelfContained    clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

set_target_properties( SelfContained    PROPERTIES VERSION ${CLPROBDIST_VERSION} )
set_target_properties( SelfContained    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

install( TARGETS SelfContained
        RUNTIME DESTINATION bin${SUFFIX_BIN}
        LIBRARY DESTINATION lib${SUFFIX_LIB}
        ARCHIVE DESTINATION lib${SUFFIX_LIB}/import
        )


if( clRNG_VERSION )

        add_executable( DocsTutorial1   ${DocsTutorial1.Files} )
        add_executable( DocsTutorial2   ${DocsTutorial2.Files} )
        add_executable( DocsTutorial3   ${DocsTutorial3.Files} )
        add_executable( DocsTutorial4   ${DocsTutorial4.Files} )
        add_executable( Inventory       ${Inventory.Files} )
        add_executable( AsianOption     ${AsianOption.Files} )
        add_executable( Performance     ${Performance.Files} )

        include_directories(   DocsTutorial1               ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( DocsTutorial1    clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        include_directories(   DocsTutorial2               ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( DocsTutorial2    clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        include_directories(   DocsTutorial3               ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( DocsTutorial3    clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        include_directories(   DocsTutorial4               ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( DocsTutorial4    clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        include_directories(   Inventory                   ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( Inventory        clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        include_directories(   AsianOption                 ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( AsianOption      clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        include_directories(   Performance                 ${clRNG_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include ../include )
        target_link_libraries( Performance      clProbDist clRNG                 ${OPENCL_LIBRARIES}    ${DL_LIB} ${MATH_LIB} )

        set_target_properties( DocsTutorial1    PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( DocsTutorial1    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        set_target_properties( DocsTutorial2    PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( DocsTutorial2    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        set_target_properties( DocsTutorial3    PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( DocsTutorial3    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        set_target_properties( DocsTutorial4    PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( DocsTutorial4    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        set_target_properties( Inventory        PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( Inventory        PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        set_target_properties( AsianOption      PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( AsianOption      PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        set_target_properties( Performance      PROPERTIES VERSION ${CLPROBDIST_VERSION} )
        set_target_properties( Performance      PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/staging" )

        install( TARGETS DocsTutorial1 DocsTutorial2 DocsTutorial3 DocsTutorial4 Inventory AsianOption Performance
                RUNTIME DESTINATION bin${SUFFIX_BIN}
                LIBRARY DESTINATION lib${SUFFIX_LIB}
                ARCHIVE DESTINATION lib${SUFFIX_LIB}/import
                )

        install( FILES
                 "AsianOption/AsianOptionKernel.cl"
                 DESTINATION
                 "./client/AsianOption" )

        install( FILES
                 "Performance/PerformanceKernels.cl"
                 DESTINATION
                 "./client/Performance" )

        install( FILES
                 "Inventory/InventoryKernels.cl"
                 DESTINATION
                 "./client/Inventory" )

else()

        message(STATUS "not building client examples that require clRNG")

endif()
