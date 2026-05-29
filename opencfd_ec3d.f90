module Flow_Var
  use precision_EC
  real(PRE_EC), save, pointer, dimension(:, :, :)::  d, uu, v, w, T, p, cc ! 密度、x-速度、y-速度、z-速度、压力、声速
  real(PRE_EC), save, pointer, dimension(:, :, :, :):: Flux                 ! i- ,j-及k-方向的通量
  real(PRE_EC), save, pointer, dimension(:, :, :):: Lvi, Lvj, Lvk, Lci, Lcj, Lck   ! 无粘项及粘性项Jocabian的谱半径

end module Flow_Var

!------------------------------------------------------------------------------------------
! 主程序 主程序 主程序
!-----------------------------------------------------------------------------------------
program main
  use Global_Var
  implicit none

  call Init_mpi

  if (my_id .eq. 0) then
    print *, "----------------- OpenCFD-EC3D ver 1.14a (MPI-OpenMP version)------------------"
    print *, "        Copyright by Li Xinliang, lixl@imech.ac.cn                            "
    print *, "        Programming by Li Xinliang  2016-1                                   "
    print *, "----------------------------------------------------------------------------- "
  end if

  call read_parameter                     ! 读取流动参数及控制信息
!$ call omp_set_num_threads(NUM_THREADS)   ! 设置OpenMP的运行线程数 （并行数目）， 本语句对openmp编译器不是注释!

!$ if (my_id == 0) then         ! 测试一下运行的进程 （openmp编译时，不是注释）
!$OMP Parallel
!$  print *, "omp run ..."
!$OMP END parallel
!$ end if

  allocate (Mesh(Num_Mesh))                                 ! 主数据结构： “网格” （其成员是“网格块”）

  if (my_id .eq. 0) call check_mesh_multigrid               ! 检查网格配置所允许的最大重数,并设定多重网格的重数

  call Init                               ! 初始化变量（分配内存，读取网格）
  call set_control_para                   ! 设定各重网格上的控制信息（数值方法、通量技术、湍流模型、时间推进方式）
  call check_mesh_quality                 ! 检查网格质量,在网格质量差的区域降低局部时间步长
  call Init_flow                          ! 初始化流场 （初值）

  if (my_id .eq. 0) print *, " Start ......"

!------------------------------------------------------------------------
! 时间推进，采用单重网格、二重网格或三重网格； 采用1阶Euler或3阶RK
  do while (Mesh(1)%tt .lt. t_end)
    call show_Wall_time()

    if (Num_Mesh .eq. 1) then                          ! 单重网格推进1个时间步
      call NS_Time_advance(1)
    else if (Num_Mesh .eq. 2) then                    ! 2重网格推进1个时间步
      call NS_2stge_multigrid
    else                                               ! 3重网格推进1个时间步
      call NS_3stge_multigrid
    end if

    !  滤波 ,可以增强稳定性. 如Kstep_Filter=0则不使用滤波
    if (Kstep_smooth .gt. 0) then
      if (mod(Mesh(1)%Kstep, Kstep_smooth) .eq. 0) call Filtering_oneMesh(1)                      ! 滤波
    end if

!  每隔一定步数输出气动力及残差（输出到屏幕及文件: force.log, Residual.dat）
    if (mod(Mesh(1)%Kstep, Kstep_show) .eq. 0) then
      call comput_force
      call output_Res(1)
    end if
! 每隔一定步数输出数据文件(flow3d.dat, PLOT3D 格式)
    if (mod(Mesh(1)%Kstep, Kstep_Save) .eq. 0) then
      call output_flow
!         if(If_debug == 1 ) call output_vt                ! 输出湍流粘性系数，供debug使用   ! Bug 2017-5-11
      if (If_debug == 1 .and. If_viscous == 1 .and. Iflag_turbulence_model .ne. 0) call output_vt                ! 输出湍流粘性系数，供debug使用
    end if

    ! 进行时间平均
    if (Kstep_average > 0) then
      if (mod(Mesh(1)%Kstep, Kstep_average) .eq. 0) then
        call Time_average           ! 时间平均
      end if
      if (mod(Mesh(1)%Kstep, Kstep_Save) .eq. 0) then
        call output_flow_average   ! 输出时均场,PLOT3D格式
      end if
    end if

  end do

end
!----------------------------------------------------------------------------------------------
subroutine Init_mpi
  use Global_var
  implicit none
  integer, parameter:: IBuffer_Size = 10000000
!      real*8,allocatable,dimension(:)::  Buffer_mpi    ! Buffer for MPI  message transfer (used by MPI_Bsend)
  real*8, save :: Buffer_mpi(IBuffer_Size)
  integer ierr

!------------------------------------------------
  call mpi_init(ierr)                                     ! 初始化MPI
  call mpi_comm_rank(MPI_COMM_WORLD, my_id, ierr)           ! 获取本进程编号
  call mpi_comm_size(MPI_COMM_WORLD, Total_proc, ierr)
!       allocate(Buffer_mpi(IBuffer_Size))
  call MPI_BUFFER_ATTACH(Buffer_mpi, 8*IBuffer_Size, ierr)   ! 创建消息发送缓冲区，供MPI_Bsend()使用
end subroutine Init_mpi

!  显示（墙钟）时间，用于统计MPI并行效率
subroutine show_Wall_time()
  use Global_var
  real*8:: wtime
  real*8, save:: wtime0, wtime1    ! 初始时间，上一步的时间
  integer, save:: KP = 0  ! 计算步
  if (my_id .eq. 0) then
    wtime = MPI_Wtime()
    if (KP .eq. 0) then
      wtime0 = wtime   ! 初始CPU时间
    else
      if (mod(Mesh(1)%Kstep, Kstep_Show) .eq. 0) then
        print *, "CPU wall time in this step:", wtime-wtime1
        print *, "Averaged CPU wall time is:", (wtime-wtime0)/KP
      end if
    end if
    wtime1 = wtime
    KP = KP+1    ! 统计计算步
  end if

end subroutine show_wall_time