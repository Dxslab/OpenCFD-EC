program add_noise
  implicit none
  integer, parameter:: PRE_EC = 8
  
  real(PRE_EC):: Ma, Re, AoA, AoS, p_outlet, t_end, gamma, PrL, PrT
  real(PRE_EC):: dt_global, CFL, dtmax, dtmin, w_LU, T_inf, Twall, Kt_inf, Wt_inf
  real(PRE_EC):: Step_Inner_Limit, Res_Inner_Limit, Res_Inner_Rel, MUT_MAX, Ref_S, Ref_L, Centroid(3)
  real(PRE_EC):: Ldmin, Ldmax, Lpmin, Lpmax, Lumax, LSAmax, CP1_NSA, CP2_NSA
  real(PRE_EC):: Turbo_P0, Turbo_T0, Turbo_L0, Turbo_w, Turbo_Periodic_seta
  real(PRE_EC):: Periodic_dX, Periodic_dY, Periodic_dZ
  integer:: Kstep_save, Iflag_savefile
  integer:: Iflag_turbulence_model, Iflag_init, If_viscous
  integer:: Iflag_local_dt, Time_Method, If_Residual_smoothing, If_dtime_mesh
  integer:: Iflag_Scheme, Iflag_Flux, IFlag_Reconstruction
  integer:: Mesh_File_Format, Kstep_show, Kstep_average, Kstep_smooth, Kstep_init_smooth
  integer:: Num_Mesh, IF_Debug, NUM_THREADS, Bound_Scheme
  integer:: Pre_Step_Mesh(10), Cood_Y_UP, IFLAG_LIMIT_FLOW, Pdebug(4)
  integer:: IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity
  integer:: IF_Innerflow
  
  namelist /control_ec/ Ma, Re, AoA, AoS, p_outlet, t_end, &
    gamma, PrL, PrT, Kstep_save, &
    Iflag_turbulence_model, Iflag_init, If_viscous, &
    Iflag_local_dt, dt_global, CFL, dtmax, dtmin, Time_Method, &
    If_Residual_smoothing, w_LU, If_dtime_mesh, &
    Iflag_Scheme, Iflag_Flux, IFlag_Reconstruction, &
    Mesh_File_Format, Kstep_show, Kstep_average, Kstep_smooth, Kstep_init_smooth, &
    Num_Mesh, T_inf, T_inf, Twall, Kt_inf, Wt_inf, IF_Debug, NUM_THREADS, &
    Step_Inner_Limit, Res_Inner_Limit, Res_Inner_Rel, MUT_MAX, Bound_Scheme, &
    Pre_Step_Mesh, Ref_S, Ref_L, Centroid, Cood_Y_UP, IFLAG_LIMIT_FLOW, Pdebug, &
    Ldmin, Ldmax, Lpmin, Lpmax, Lumax, LSAmax, CP1_NSA, CP2_NSA, &
    IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, &
    Turbo_P0, Turbo_T0, Turbo_L0, Turbo_w, Turbo_Periodic_seta, &
    Periodic_dX, Periodic_dY, Periodic_dZ, &
    IF_Innerflow, Iflag_savefile
    
  integer:: Num_Block, k, m, m1, nx, ny, nz, i, j
  integer, allocatable:: NI(:), NJ(:), NK(:)
  real(PRE_EC), allocatable:: U(:,:,:,:)
  real(PRE_EC):: rand_val
  integer:: tg_block, i1, i2, j1, j2, k1, k2
  real(PRE_EC):: du, dv, dw, mag
  integer:: NVAR1

  !---- default--------
  Mesh_File_Format = 0
  Iflag_turbulence_model = 0
  
  open (99, file="control.ec", status="old")
  read (99, nml=control_ec)
  close (99)
  
  if (Iflag_turbulence_model==0 .or. Iflag_turbulence_model==1) then
     NVAR1 = 5
  else if (Iflag_turbulence_model==2) then
     NVAR1 = 6
  else if (Iflag_turbulence_model==3) then
     NVAR1 = 7
  end if

  if (Mesh_File_Format .eq. 0) then
    open (99, file="Mesh3d.dat", form="unformatted")
    read (99) Num_Block
  else
    open (99, file="Mesh3d.dat")
    read (99, *) Num_Block
  end if
  
  allocate (NI(Num_Block), NJ(Num_Block), NK(Num_Block))
  
  if (Mesh_File_Format .eq. 0) then
    read (99) (NI(k), NJ(k), NK(k), k=1, Num_Block)
  else
    read (99, *) (NI(k), NJ(k), NK(k), k=1, Num_Block)
  end if
  close (99)

  print *, "Number of Blocks:", Num_Block
  
  print *, "Enter target block ID to perturb (1 to ", Num_Block, "):"
  read *, tg_block
  nx = NI(tg_block)
  ny = NJ(tg_block)
  nz = NK(tg_block)
  print *, "Block dimensions: nx=", nx, " ny=", ny, " nz=", nz
  
  print *, "Enter index range to perturb: i1, i2, j1, j2, k1, k2"
  print *, "e.g., 5 10 5 10 5 10 (Max limits are ", nx, ny, nz, ")"
  read *, i1, i2, j1, j2, k1, k2
  
  print *, "Enter max perturbation magnitude for U, V, W (e.g., 0.1 0.1 0.1):"
  read *, du, dv, dw

  open (99, file="flow3d.dat", form="unformatted", status="old")
  open (100, file="flow3d_new.dat", form="unformatted")
  
  do m = 1, Num_Block
    nx = NI(m)
    ny = NJ(m)
    nz = NK(m)
    allocate(U(0:nx, 0:ny, 0:nz, 1:5))
    read (99) ((((U(i, j, k, m1), i=0, nx), j=0, ny), k=0, nz), m1=1, 5)
    
    if (m == tg_block) then
      print *, "Applying perturbations to Block ", m
      do k = max(0,k1), min(nz,k2)
        do j = max(0,j1), min(ny,j2)
          do i = max(0,i1), min(nx,i2)
            call random_number(rand_val)
            ! asymmetric random between 0 and 1 (so mean is non-zero, adding asymmetry)
            U(i,j,k,2) = U(i,j,k,2) + du * rand_val
            call random_number(rand_val)
            U(i,j,k,3) = U(i,j,k,3) + dv * rand_val
            call random_number(rand_val)
            U(i,j,k,4) = U(i,j,k,4) + dw * rand_val
          end do
        end do
      end do
    end if
    
    write (100) ((((U(i, j, k, m1), i=0, nx), j=0, ny), k=0, nz), m1=1, 5)
    deallocate(U)
  end do
  
  close(99)
  close(100)
  
  ! Also copy SA3d / SST3d if they exist
  if (NVAR1 == 6) then
     open (99, file="SA3d.dat", form="unformatted", status="old", iostat=k)
     if (k==0) then
       open(100, file="SA3d_new.dat", form="unformatted")
       do m = 1, Num_Block
         nx = NI(m); ny = NJ(m); nz = NK(m)
         allocate(U(0:nx, 0:ny, 0:nz, 1:1))
         read (99) (((U(i, j, k, 1), i=0, nx), j=0, ny), k=0, nz)
         write (100) (((U(i, j, k, 1), i=0, nx), j=0, ny), k=0, nz)
         deallocate(U)
       end do
       close(99)
       close(100)
       print *, "Processed SA3d.dat"
     end if
  end if
  if (NVAR1 == 7) then
     open (99, file="SST3d.dat", form="unformatted", status="old", iostat=k)
     if (k==0) then
       open(100, file="SST3d_new.dat", form="unformatted")
       do m = 1, Num_Block
         nx = NI(m); ny = NJ(m); nz = NK(m)
         allocate(U(0:nx, 0:ny, 0:nz, 1:2))
         read (99) ((((U(i, j, k, m1), i=0, nx), j=0, ny), k=0, nz), m1=1, 2)
         write (100) ((((U(i, j, k, m1), i=0, nx), j=0, ny), k=0, nz), m1=1, 2)
         deallocate(U)
       end do
       close(99)
       close(100)
       print *, "Processed SST3d.dat"
     end if
  end if
  
  print *, "Done. New flow files generated as *_new.dat"

end program add_noise
