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
add_executable(readflow3d-tecio util/readflow3d-tecio.f90)

# Link MPI and OpenMP to the executable
target_link_libraries(opencfd-ec PRIVATE MPI::MPI_Fortran OpenMP::OpenMP_Fortran)
target_include_directories(readflow3d-tecio PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/util/libtecio)
target_link_libraries(readflow3d-tecio PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/util/libtecio/libtecio.a
    pthread
    stdc++
)

# Enable LTO (Interprocedural Optimization) which corresponds to -flto or -ipo
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if(ipo_supported)
    set_property(TARGET opencfd-ec PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET readflow3d-tecio PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Compiler specific optimization and flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # For gfortran
    target_compile_options(opencfd-ec PRIVATE
        -O3
        -march=native
        -ffree-line-length-none # Required because some lines exceed 132 columns
    )
    target_compile_options(readflow3d-tecio PRIVATE
        -O3
        -march=native
        -fcray-pointer
        -ffree-line-length-none # Required because some lines exceed 132 columns
    )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    # For ifort or ifx
    target_compile_options(opencfd-ec PRIVATE
        -O3
        -xHost
        -fp-model precise
    )
    target_compile_options(readflow3d-tecio PRIVATE
        -O3
        -xHost
        -fp-model precise
    )
endif()
