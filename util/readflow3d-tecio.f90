!-----------------------------------------------------------
! Direct Tecplot binary exporter based on readflow3d-ver2.5
! Reads Mesh3d.dat + flow3d.dat/flow3d_average.dat and writes
! the full 3D flow field through libtecio.a.
!-----------------------------------------------------------
module Tecio_Flow_Data
  implicit none

  integer, parameter :: PRE_EC = 8

  real(PRE_EC), save :: Ma = 1.0d0
  real(PRE_EC), save :: Re = 100.0d0
  real(PRE_EC), save :: gamma = 1.4d0
  real(PRE_EC), save :: T_inf = 288.15d0
  real(PRE_EC), save :: p00 = 0.0d0
  integer, save :: Num_Block = 0
  integer, save :: Mesh_File_Format = 0
  character(len=256), save :: Flow_File = 'flow3d.dat'
  character(len=256), save :: Tecplot_File = 'flow3d-tecio.plt'

  type Block_TYPE
    integer :: nx = 0
    integer :: ny = 0
    integer :: nz = 0
    real(PRE_EC), allocatable :: xc(:, :, :), yc(:, :, :), zc(:, :, :)
    real(PRE_EC), allocatable :: dc(:, :, :), uc(:, :, :), vc(:, :, :), wc(:, :, :)
    real(PRE_EC), allocatable :: Tc(:, :, :), pc(:, :, :)
    real(PRE_EC), allocatable :: mut(:, :, :), vt(:, :, :), Kt(:, :, :), Wt(:, :, :)
    real(PRE_EC), allocatable :: x(:, :, :), y(:, :, :), z(:, :, :)
    real(PRE_EC), allocatable :: d(:, :, :), u(:, :, :), v(:, :, :), w(:, :, :)
    real(PRE_EC), allocatable :: T(:, :, :), p(:, :, :), mut1(:, :, :), vt1(:, :, :)
    real(PRE_EC), allocatable :: Kt1(:, :, :), Wt1(:, :, :)
    real(PRE_EC), allocatable :: dw(:, :, :)
  end type Block_TYPE

  type(Block_TYPE), allocatable, target, save :: Block(:)
end module Tecio_Flow_Data

program readflow3d_tecio
  use Tecio_Flow_Data
  implicit none

  include 'tecio.f90'

  call parse_command_line()
  call init_flow_data()
  call write_tecplot_binary()

contains

  subroutine parse_command_line()
    implicit none

    integer :: nargs
    character(len=256) :: arg

    nargs = command_argument_count()
    if (nargs >= 1) then
      call get_command_argument(1, arg)
      if (trim(arg) == '-h' .or. trim(arg) == '--help') then
        call print_usage()
        stop
      end if
      Flow_File = trim(arg)
    end if

    if (nargs >= 2) then
      call get_command_argument(2, Tecplot_File)
    end if
  end subroutine parse_command_line

  subroutine print_usage()
    implicit none

    print *, 'Usage: readflow3d-tecio [flow_file] [output_plt]'
    print *, 'Default flow_file : flow3d.dat'
    print *, 'Default output    : flow3d-tecio.plt'
    print *, 'Example           : ./readflow3d-tecio flow3d_average.dat flow3d_average.plt'
  end subroutine print_usage

  subroutine init_flow_data()
    implicit none

    call read_parameter_ec()
    p00 = 1.0d0/(gamma*Ma*Ma)

    print *, 'Ma, Re, gamma, T_inf = ', Ma, Re, gamma, T_inf
    print *, 'Mesh_File_Format     = ', Mesh_File_Format
    print *, 'Flow field file      = ', trim(Flow_File)
    print *, 'Tecplot output file  = ', trim(Tecplot_File)

    call read_mesh_file()
    call read_wall_distance()
    call read_flow_file(trim(Flow_File))
    call read_turbulence_files()
    call compute_value_in_mesh()
  end subroutine init_flow_data

  subroutine read_mesh_file()
    implicit none

    integer :: i, j, k, m
    integer :: nx, ny, nz
    integer, allocatable :: NI(:), NJ(:), NK(:)
    type(Block_TYPE), pointer :: B

    print *, 'Read Mesh3d.dat ...'
    if (Mesh_File_Format == 0) then
      open (99, file='Mesh3d.dat', form='unformatted', status='old', action='read')
      read (99) Num_Block
    else
      open (99, file='Mesh3d.dat', status='old', action='read')
      read (99, *) Num_Block
    end if

    print *, 'Num_Block = ', Num_Block

    allocate (Block(Num_Block))
    allocate (NI(Num_Block), NJ(Num_Block), NK(Num_Block))

    if (Mesh_File_Format == 0) then
      read (99) (NI(k), NJ(k), NK(k), k=1, Num_Block)
    else
      read (99, *) (NI(k), NJ(k), NK(k), k=1, Num_Block)
    end if

    do m = 1, Num_Block
      B => Block(m)
      B%nx = NI(m)
      B%ny = NJ(m)
      B%nz = NK(m)
      nx = B%nx
      ny = B%ny
      nz = B%nz

      allocate (B%xc(0:nx, 0:ny, 0:nz), B%yc(0:nx, 0:ny, 0:nz), B%zc(0:nx, 0:ny, 0:nz))
      allocate (B%dc(0:nx, 0:ny, 0:nz), B%uc(0:nx, 0:ny, 0:nz), B%vc(0:nx, 0:ny, 0:nz), &
                B%wc(0:nx, 0:ny, 0:nz), B%Tc(0:nx, 0:ny, 0:nz), B%pc(0:nx, 0:ny, 0:nz), &
                B%mut(0:nx, 0:ny, 0:nz), B%vt(0:nx, 0:ny, 0:nz), B%Kt(0:nx, 0:ny, 0:nz), &
                B%Wt(0:nx, 0:ny, 0:nz))
      allocate (B%x(nx, ny, nz), B%y(nx, ny, nz), B%z(nx, ny, nz))
      allocate (B%d(nx, ny, nz), B%u(nx, ny, nz), B%v(nx, ny, nz), B%w(nx, ny, nz), &
                B%T(nx, ny, nz), B%p(nx, ny, nz), B%mut1(nx, ny, nz), B%vt1(nx, ny, nz), &
                B%Kt1(nx, ny, nz), B%Wt1(nx, ny, nz))
      allocate (B%dw(nx-1, ny-1, nz-1))

      B%mut = 0.0d0
      B%vt = 0.0d0
      B%Kt = 0.0d0
      B%Wt = 0.0d0
      B%mut1 = 0.0d0
      B%vt1 = 0.0d0
      B%Kt1 = 0.0d0
      B%Wt1 = 0.0d0
      B%dw = 0.0d0

      print *, '  Block ', m, ' size = ', nx, ny, nz
      if (Mesh_File_Format == 0) then
        read (99) (((B%x(i, j, k), i=1, nx), j=1, ny), k=1, nz), &
                  (((B%y(i, j, k), i=1, nx), j=1, ny), k=1, nz), &
                  (((B%z(i, j, k), i=1, nx), j=1, ny), k=1, nz)
      else
        read (99, *) (((B%x(i, j, k), i=1, nx), j=1, ny), k=1, nz), &
                     (((B%y(i, j, k), i=1, nx), j=1, ny), k=1, nz), &
                     (((B%z(i, j, k), i=1, nx), j=1, ny), k=1, nz)
      end if

      do k = 1, nz-1
        do j = 1, ny-1
          do i = 1, nx-1
            B%xc(i, j, k) = (B%x(i, j, k)+B%x(i+1, j, k)+B%x(i, j+1, k)+ &
                             B%x(i, j, k+1)+B%x(i+1, j+1, k)+B%x(i+1, j, k+1)+ &
                             B%x(i, j+1, k+1)+B%x(i+1, j+1, k+1))/8.0d0
            B%yc(i, j, k) = (B%y(i, j, k)+B%y(i+1, j, k)+B%y(i, j+1, k)+ &
                             B%y(i, j, k+1)+B%y(i+1, j+1, k)+B%y(i+1, j, k+1)+ &
                             B%y(i, j+1, k+1)+B%y(i+1, j+1, k+1))/8.0d0
            B%zc(i, j, k) = (B%z(i, j, k)+B%z(i+1, j, k)+B%z(i, j+1, k)+ &
                             B%z(i, j, k+1)+B%z(i+1, j+1, k)+B%z(i+1, j, k+1)+ &
                             B%z(i, j+1, k+1)+B%z(i+1, j+1, k+1))/8.0d0
          end do
        end do
      end do
    end do

    close (99)
  end subroutine read_mesh_file

  subroutine read_wall_distance()
    implicit none

    integer :: i, j, k, m
    logical :: exists
    type(Block_TYPE), pointer :: B

    inquire (file='wall_dist.dat', exist=exists)
    if (.not. exists) then
      print *, 'wall_dist.dat not found, skip.'
      return
    end if

    print *, 'Read wall_dist.dat ...'
    open (100, file='wall_dist.dat', form='unformatted', status='old', action='read')
    do m = 1, Num_Block
      B => Block(m)
      if (B%nx > 1 .and. B%ny > 1 .and. B%nz > 1) then
        read (100) (((B%dw(i, j, k), i=1, B%nx-1), j=1, B%ny-1), k=1, B%nz-1)
      end if
    end do
    close (100)
  end subroutine read_wall_distance

  subroutine read_flow_file(filename)
    implicit none

    character(len=*), intent(in) :: filename
    integer :: i, j, k, m
    type(Block_TYPE), pointer :: B

    open (99, file=trim(filename), form='unformatted', status='old', action='read')
    do m = 1, Num_Block
      B => Block(m)
      read (99) (((B%dc(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz), &
                (((B%uc(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz), &
                (((B%vc(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz), &
                (((B%wc(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz), &
                (((B%Tc(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz)

      do k = 0, B%nz
        do j = 0, B%ny
          do i = 0, B%nx
            B%pc(i, j, k) = p00*B%dc(i, j, k)*B%Tc(i, j, k)
          end do
        end do
      end do
    end do
    close (99)
  end subroutine read_flow_file

  subroutine read_turbulence_files()
    implicit none

    integer :: i, j, k, m
    logical :: exists
    type(Block_TYPE), pointer :: B

    inquire (file='SA3d.dat', exist=exists)
    if (exists) then
      print *, 'Read SA3d.dat ...'
      open (100, file='SA3d.dat', form='unformatted', status='old', action='read')
      do m = 1, Num_Block
        B => Block(m)
        read (100) (((B%vt(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz)
      end do
      close (100)
    else
      print *, 'SA3d.dat not found, vt set to zero.'
    end if

    inquire (file='SST3d.dat', exist=exists)
    if (exists) then
      print *, 'Read SST3d.dat ...'
      open (100, file='SST3d.dat', form='unformatted', status='old', action='read')
      do m = 1, Num_Block
        B => Block(m)
        read (100) (((B%Kt(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz), &
                   (((B%Wt(i, j, k), i=0, B%nx), j=0, B%ny), k=0, B%nz)
      end do
      close (100)
    else
      print *, 'SST3d.dat not found, Kt/Wt set to zero.'
    end if
  end subroutine read_turbulence_files

  subroutine compute_value_in_mesh()
    implicit none

    integer :: i, j, k, m
    type(Block_TYPE), pointer :: B

    do m = 1, Num_Block
      B => Block(m)
      do k = 1, B%nz
        do j = 1, B%ny
          do i = 1, B%nx
            B%d(i, j, k) = (B%dc(i-1, j-1, k-1)+B%dc(i, j-1, k-1)+B%dc(i-1, j, k-1)+B%dc(i, j, k-1) &
                            +B%dc(i-1, j-1, k)+B%dc(i, j-1, k)+B%dc(i-1, j, k)+B%dc(i, j, k))/8.0d0
            B%u(i, j, k) = (B%uc(i-1, j-1, k-1)+B%uc(i, j-1, k-1)+B%uc(i-1, j, k-1)+B%uc(i, j, k-1) &
                            +B%uc(i-1, j-1, k)+B%uc(i, j-1, k)+B%uc(i-1, j, k)+B%uc(i, j, k))/8.0d0
            B%v(i, j, k) = (B%vc(i-1, j-1, k-1)+B%vc(i, j-1, k-1)+B%vc(i-1, j, k-1)+B%vc(i, j, k-1) &
                            +B%vc(i-1, j-1, k)+B%vc(i, j-1, k)+B%vc(i-1, j, k)+B%vc(i, j, k))/8.0d0
            B%w(i, j, k) = (B%wc(i-1, j-1, k-1)+B%wc(i, j-1, k-1)+B%wc(i-1, j, k-1)+B%wc(i, j, k-1) &
                            +B%wc(i-1, j-1, k)+B%wc(i, j-1, k)+B%wc(i-1, j, k)+B%wc(i, j, k))/8.0d0
            B%T(i, j, k) = (B%Tc(i-1, j-1, k-1)+B%Tc(i, j-1, k-1)+B%Tc(i-1, j, k-1)+B%Tc(i, j, k-1) &
                            +B%Tc(i-1, j-1, k)+B%Tc(i, j-1, k)+B%Tc(i-1, j, k)+B%Tc(i, j, k))/8.0d0
            B%p(i, j, k) = (B%pc(i-1, j-1, k-1)+B%pc(i, j-1, k-1)+B%pc(i-1, j, k-1)+B%pc(i, j, k-1) &
                            +B%pc(i-1, j-1, k)+B%pc(i, j-1, k)+B%pc(i-1, j, k)+B%pc(i, j, k))/8.0d0
            B%mut1(i, j, k) = (B%mut(i-1, j-1, k-1)+B%mut(i, j-1, k-1)+B%mut(i-1, j, k-1)+B%mut(i, j, k-1) &
                               +B%mut(i-1, j-1, k)+B%mut(i, j-1, k)+B%mut(i-1, j, k)+B%mut(i, j, k))/8.0d0
            B%vt1(i, j, k) = (B%vt(i-1, j-1, k-1)+B%vt(i, j-1, k-1)+B%vt(i-1, j, k-1)+B%vt(i, j, k-1) &
                              +B%vt(i-1, j-1, k)+B%vt(i, j-1, k)+B%vt(i-1, j, k)+B%vt(i, j, k))/8.0d0
            B%Kt1(i, j, k) = (B%Kt(i-1, j-1, k-1)+B%Kt(i, j-1, k-1)+B%Kt(i-1, j, k-1)+B%Kt(i, j, k-1) &
                              +B%Kt(i-1, j-1, k)+B%Kt(i, j-1, k)+B%Kt(i-1, j, k)+B%Kt(i, j, k))/8.0d0
            B%Wt1(i, j, k) = (B%Wt(i-1, j-1, k-1)+B%Wt(i, j-1, k-1)+B%Wt(i-1, j, k-1)+B%Wt(i, j, k-1) &
                              +B%Wt(i-1, j-1, k)+B%Wt(i, j-1, k)+B%Wt(i-1, j, k)+B%Wt(i, j, k))/8.0d0
          end do
        end do
      end do
    end do
  end subroutine compute_value_in_mesh

  subroutine write_tecplot_binary()
    implicit none

    integer :: m, npts
    real(kind=8), allocatable :: mut_out(:, :, :)
    character(len=128) :: dataset_title
    character(len=64) :: var_list
    character(len=64) :: zone_title
    character(len=256) :: output_name
    character*1 :: NULCHAR
    integer*4 :: Debug, FileType, FileFormat, VIsDouble
    integer*4 :: I_ret, IMax, JMax, KMax, ZoneType
    integer*4 :: ICellMax, JCellMax, KCellMax, IsBlock
    integer*4 :: NFConns, FNMode, TotalNumFaceNodes
    integer*4 :: TotalNumBndryFaces, TotalNumBndryConnections, ShrConn
    integer*4 :: StrandID, ParentZn, III
    integer*8 :: NullPtr
    integer*4 :: Null(*)
    real*8 :: SolTime
    pointer (NullPtr, Null)

    dataset_title = 'OpenCFD-EC flow field'
    var_list = 'X Y Z Density U V W Temperature Pressure mut vt Kt Wt'
    output_name = trim(Tecplot_File)

    Debug = 1
    VIsDouble = 1
    FileType = 0
    FileFormat = 0
    NULCHAR = CHAR(0)
    NullPtr = 0

    I_ret = TECINI142(trim(dataset_title)//NULCHAR, &
                      trim(var_list)//NULCHAR, &
                      trim(output_name)//NULCHAR, &
                      '.'//NULCHAR, &
                      FileFormat, FileType, Debug, VIsDouble)
    call check_tecio_status(I_ret, 'TECINI142')

    do m = 1, Num_Block
      IMax = Block(m)%nx
      JMax = Block(m)%ny
      KMax = Block(m)%nz
      ZoneType = 0
      SolTime = 0.0d0
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      TotalNumFaceNodes = 0
      TotalNumBndryFaces = 0
      TotalNumBndryConnections = 0
      ShrConn = 0

      write (zone_title, '(A,I0)') 'Block-', m
      I_ret = TECZNE142(trim(zone_title)//NULCHAR, &
                        ZoneType, IMax, JMax, KMax, &
                        ICellMax, JCellMax, KCellMax, &
                        SolTime, StrandID, ParentZn, IsBlock, &
                        NFConns, FNMode, &
                        TotalNumFaceNodes, TotalNumBndryFaces, TotalNumBndryConnections, &
                        Null, Null, Null, ShrConn)
      call check_tecio_status(I_ret, 'TECZNE142 for '//trim(zone_title))

      III = IMax*JMax*KMax
      allocate (mut_out(IMax, JMax, KMax))
      mut_out = Block(m)%mut1*Re

      I_ret = TECDATD142(III, Block(m)%x)
      call check_tecio_status(I_ret, 'TECDATD142 X for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%y)
      call check_tecio_status(I_ret, 'TECDATD142 Y for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%z)
      call check_tecio_status(I_ret, 'TECDATD142 Z for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%d)
      call check_tecio_status(I_ret, 'TECDATD142 Density for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%u)
      call check_tecio_status(I_ret, 'TECDATD142 U for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%v)
      call check_tecio_status(I_ret, 'TECDATD142 V for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%w)
      call check_tecio_status(I_ret, 'TECDATD142 W for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%T)
      call check_tecio_status(I_ret, 'TECDATD142 Temperature for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%p)
      call check_tecio_status(I_ret, 'TECDATD142 Pressure for '//trim(zone_title))
      I_ret = TECDATD142(III, mut_out)
      call check_tecio_status(I_ret, 'TECDATD142 mut for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%vt1)
      call check_tecio_status(I_ret, 'TECDATD142 vt for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%Kt1)
      call check_tecio_status(I_ret, 'TECDATD142 Kt for '//trim(zone_title))
      I_ret = TECDATD142(III, Block(m)%Wt1)
      call check_tecio_status(I_ret, 'TECDATD142 Wt for '//trim(zone_title))

      deallocate (mut_out)
      npts = III
      print *, '  Wrote ', trim(zone_title), ', points = ', npts
    end do

    I_ret = TECEND142()
    call check_tecio_status(I_ret, 'TECEND142')

    print *, 'Tecplot binary file written: ', trim(Tecplot_File)
  end subroutine write_tecplot_binary

  subroutine check_tecio_status(ierr, where_msg)
    implicit none

    integer*4, intent(in) :: ierr
    character(len=*), intent(in) :: where_msg

    if (ierr /= 0) then
      print *, 'TecIO call failed at ', trim(where_msg), ', ierr = ', ierr
      stop 1
    end if
  end subroutine check_tecio_status

  subroutine read_parameter_ec()
    implicit none

    real(PRE_EC), parameter :: PI = 3.14159265358979d0
    real(PRE_EC) :: R0, a0, d0, mu0, mu1
    real(PRE_EC) :: AoA, AoS, p_outlet, t_end
    real(PRE_EC) :: PrL, PrT
    real(PRE_EC) :: dt_global, CFL, dtmax, dtmin, Time_Method
    real(PRE_EC) :: w_LU, Twall, Kt_inf, Wt_inf
    real(PRE_EC) :: Step_Inner_Limit, Res_Inner_Limit, MUT_MAX
    real(PRE_EC) :: Ref_S, Ref_L, Centroid(3)
    real(PRE_EC) :: Ldmin, Ldmax, Lpmin, Lpmax, Lumax, LSAmax, CP1_NSA, CP2_NSA
    real(PRE_EC) :: Turbo_P0, Turbo_T0, Turbo_L0, Turbo_w, Turbo_Periodic_seta
    real(PRE_EC) :: Periodic_dX, Periodic_dY, Periodic_dZ
    integer :: Kstep_save, Kstep_average, Iflag_turbulence_model, Iflag_init, If_viscous
    integer :: Iflag_local_dt, If_Residual_smoothing, If_dtime_mesh
    integer :: Iflag_Scheme, Iflag_Flux, IFlag_Reconstruction, Pre_Step_Mesh(3)
    integer :: Kstep_show, Kstep_smooth, Kstep_init_smooth, Bound_Scheme
    integer :: Num_Mesh, IF_Debug, NUM_THREADS, Cood_Y_UP, IFLAG_LIMIT_FLOW, Pdebug(4)
    integer :: IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, IF_Innerflow
    integer :: Iflag_savefile

    namelist /control_ec/ Ma, Re, AoA, AoS, p_outlet, t_end, &
      gamma, PrL, PrT, &
      Kstep_save, &
      Iflag_turbulence_model, Iflag_init, If_viscous, &
      Iflag_local_dt, dt_global, CFL, dtmax, dtmin, Time_Method, &
      If_Residual_smoothing, w_LU, If_dtime_mesh, &
      Iflag_Scheme, Iflag_Flux, IFlag_Reconstruction, &
      Mesh_File_Format, Kstep_show, Kstep_average, Kstep_smooth, Kstep_init_smooth, &
      Num_Mesh, T_inf, Twall, Kt_inf, Wt_inf, IF_Debug, NUM_THREADS, &
      Step_Inner_Limit, Res_Inner_Limit, MUT_MAX, Bound_Scheme, &
      Pre_Step_Mesh, Ref_S, Ref_L, Centroid, Cood_Y_UP, IFLAG_LIMIT_FLOW, Pdebug, &
      Ldmin, Ldmax, Lpmin, Lpmax, Lumax, LSAmax, CP1_NSA, CP2_NSA, &
      IF_TurboMachinary, Ref_medium_usrdef, IF_Scheme_Positivity, &
      Turbo_P0, Turbo_T0, Turbo_L0, Turbo_w, Turbo_Periodic_seta, &
      Periodic_dX, Periodic_dY, Periodic_dZ, &
      IF_Innerflow, Iflag_savefile

    Ma = 1.0d0
    Re = 100.0d0
    gamma = 1.4d0
    T_inf = 288.15d0
    Mesh_File_Format = 0

    AoA = 0.0d0
    AoS = 0.0d0
    p_outlet = 0.0d0
    t_end = 0.0d0
    PrL = 0.7d0
    PrT = 0.9d0
    dt_global = 0.0d0
    CFL = 0.0d0
    dtmax = 0.0d0
    dtmin = 0.0d0
    Time_Method = 0.0d0
    w_LU = 0.0d0
    Twall = 0.0d0
    Kt_inf = 0.0d0
    Wt_inf = 0.0d0
    Step_Inner_Limit = 0.0d0
    Res_Inner_Limit = 0.0d0
    MUT_MAX = 0.0d0
    Ref_S = 0.0d0
    Ref_L = 0.0d0
    Centroid = 0.0d0
    Ldmin = 0.0d0
    Ldmax = 0.0d0
    Lpmin = 0.0d0
    Lpmax = 0.0d0
    Lumax = 0.0d0
    LSAmax = 0.0d0
    CP1_NSA = 0.0d0
    CP2_NSA = 0.0d0
    Turbo_P0 = 0.0d0
    Turbo_T0 = T_inf
    Turbo_L0 = 1.0d0
    Turbo_w = 0.0d0
    Turbo_Periodic_seta = 0.0d0
    Periodic_dX = 0.0d0
    Periodic_dY = 0.0d0
    Periodic_dZ = 0.0d0

    Kstep_save = 0
    Kstep_average = 0
    Iflag_turbulence_model = 0
    Iflag_init = 0
    If_viscous = 0
    Iflag_local_dt = 0
    If_Residual_smoothing = 0
    If_dtime_mesh = 0
    Iflag_Scheme = 0
    Iflag_Flux = 0
    IFlag_Reconstruction = 0
    Pre_Step_Mesh = 0
    Kstep_show = 0
    Kstep_smooth = 0
    Kstep_init_smooth = 0
    Bound_Scheme = 0
    Num_Mesh = 0
    IF_Debug = 0
    NUM_THREADS = 0
    Cood_Y_UP = 0
    IFLAG_LIMIT_FLOW = 0
    Pdebug = 0
    IF_TurboMachinary = 0
    Ref_medium_usrdef = 1
    IF_Scheme_Positivity = 0
    IF_Innerflow = 0
    Iflag_savefile = 0

    open (99, file='control.ec', status='old', action='read')
    read (99, nml=control_ec)
    close (99)

    if ((IF_TurboMachinary == 1 .or. IF_Innerflow == 1) .and. Ref_medium_usrdef == 0) then
      T_inf = Turbo_T0
      gamma = 1.4d0
      PrL = 0.7d0
      PrT = 0.9d0
      R0 = 287.06d0
      a0 = sqrt(gamma*R0*Turbo_T0)
      mu0 = 1.179d-5
      mu1 = mu0*sqrt((Turbo_T0/288.15d0)**3)*(288.15d0+110.4d0)/(Turbo_T0+110.4d0)
      d0 = Turbo_P0/(R0*Turbo_T0)
      Re = d0*a0*Turbo_L0/mu1
      Ma = 1.0d0
      Turbo_w = 2.0d0*PI*Turbo_w/(a0/Turbo_L0)
      p_outlet = p_outlet/Turbo_P0
    end if
  end subroutine read_parameter_ec

end program readflow3d_tecio
