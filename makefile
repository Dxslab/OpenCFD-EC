FC=mpifort
opt=-O3 -flto=auto -march=native -ffree-line-length-none -fopenmp # for gnu compiler
# opt=-O3 -xHost -fp-model precise -ipo -qopenmp # for intel compiler and cpu

srs=  sub_modules.f90 \
      opencfd_ec3d.f90 \
      sub_interfaces.f90 \
      sub_init.f90 \
      sub_Residual.f90 \
      sub_flux_split.f90 \
      sub_boundary.f90 \
      sub_turbulence_BL.f90 \
      sub_turbulence_SA.f90 \
      sub_turbulence_NewSA.f90 \
      sub_turbulence_SST.f90 \
      sub_time_advance.f90 \
      sub_update_buffer_mpi.f90 \
      sub_Post.f90  \
      sub_time_acceleraction.f90 \
      sub_LU_SGS.f90  \
      sub_IO.f90 \
      sub_geometry.f90 \
      sub_partation_mpi.f90 \
      sub_convert_inp.f90 \
      sub_read_parameter.f90 \
      sub_scheme.f90 \
      sub_comput_dw.f90 \
      sub_debug.f90 \
      sub_filtering.f90 \
      sub_boundary_user.f90 \
      sub_update_Meshlink_mpi.f90 \
      sub_Finite_Difference1.f90 \
      sub_Finite_Difference2.f90 \
      sub_Post_timeAverage.f90 \
      sub_limitflow.f90


OBJS=$(srs:.f90=.o)

%.o:%.f90
	$(FC) $(opt) -c $<


default: $(OBJS)
	$(FC) $(opt)  -o opencfd-ec.out $(OBJS)

clean:
	rm -f *.out *.o work.* *.pc *.mod

