cmake_minimum_required(VERSION 3.10)
project(OpenCFD_EC LANGUAGES Fortran)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(SOURCES
    sub_modules.f90
    opencfd_ec3d.f90
    sub_interfaces.f90
    sub_init.f90
    sub_Residual.f90
    sub_flux_split.f90
    sub_boundary.f90
    sub_turbulence_BL.f90
    sub_turbulence_SA.f90
    sub_turbulence_NewSA.f90
    sub_turbulence_SST.f90
    sub_time_advance.f90
    sub_update_buffer_mpi.f90
    sub_Post.f90
    sub_time_acceleraction.f90
    sub_LU_SGS.f90
    sub_IO.f90
    sub_geometry.f90
    sub_partation_mpi.f90
    sub_convert_inp.f90
    sub_read_parameter.f90
    sub_scheme.f90
    sub_comput_dw.f90
    sub_debug.f90
    sub_filtering.f90
    sub_boundary_user.f90
    sub_update_Meshlink_mpi.f90
    sub_Finite_Difference1.f90
    sub_Finite_Difference2.f90
    sub_Post_timeAverage.f90
    sub_limitflow.f90
)

add_executable(opencfd-ec ${SOURCES})
# Keep solver .mod files out of the build root: gfortran searches the current
# directory first, so util tools would pick up same-named solver modules (e.g. Type_def1)
set_target_properties(opencfd-ec PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/opencfd-ec)

# Standalone pre/post-processing tools in util/ (one executable per source file)
set(UTIL_TOOLS
    Block_cut_1.1
    In2inp-1.1
    PlotZplane-1.3a
    add_noise
    check-mesh-ver1.1
    comput-Q-ec1.0
    convertMesh2d-3d-v1.1
    convertXYZ-ver1.1
    get_bcin_v1.4.2
    partation-1.2
    readflow3d-ver2.5
    readflow3d-wing
    readflow3d-tecio
)

# convert-cgns-inp needs the CGNS library; build it only when CGNS is found
find_path(CGNS_INCLUDE_DIR cgnslib_f.h)
find_library(CGNS_LIBRARY cgns)
if(CGNS_INCLUDE_DIR AND CGNS_LIBRARY)
    list(APPEND UTIL_TOOLS convert-cgns-inp-1.2)
else()
    message(STATUS "CGNS not found, skipping util/convert-cgns-inp-1.2.f90")
endif()

foreach(tool ${UTIL_TOOLS})
    add_executable(${tool} util/${tool}.f90)
    # Several tools define modules with the same name, keep their .mod files apart
    set_target_properties(${tool} PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/${tool})
endforeach()

# Link MPI and OpenMP to the executable
target_link_libraries(opencfd-ec PRIVATE MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
target_include_directories(readflow3d-tecio PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/util/libtecio)
target_link_libraries(readflow3d-tecio PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/util/libtecio/libtecio.a
    pthread
    stdc++
)
if(TARGET convert-cgns-inp-1.2)
    target_include_directories(convert-cgns-inp-1.2 PRIVATE ${CGNS_INCLUDE_DIR})
    target_link_libraries(convert-cgns-inp-1.2 PRIVATE ${CGNS_LIBRARY})
endif()

# Enable LTO (Interprocedural Optimization) which corresponds to -flto or -ipo
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if(ipo_supported)
    set_property(TARGET opencfd-ec PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET ${UTIL_TOOLS} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Compiler specific optimization and flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # For gfortran
    target_compile_options(opencfd-ec PRIVATE
        -O3
        -march=native
        -ffree-line-length-none # Required because some lines exceed 132 columns
    )
    foreach(tool ${UTIL_TOOLS})
        target_compile_options(${tool} PRIVATE
            -O3
            -march=native
            -ffree-line-length-none
        )
    endforeach()
    target_compile_options(readflow3d-tecio PRIVATE -fcray-pointer)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    # For ifort or ifx
    target_compile_options(opencfd-ec PRIVATE
        -O3
        -xHost
        -fp-model precise
    )
    foreach(tool ${UTIL_TOOLS})
        target_compile_options(${tool} PRIVATE
            -O3
            -xHost
            -fp-model precise
        )
    endforeach()
endif()
