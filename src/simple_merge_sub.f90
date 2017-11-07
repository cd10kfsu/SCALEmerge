program main
!-------------------------------------------------------------------------------
! same as simple_merge.f90 except that it uses the subroutines to get
! the whole domain field, instead of combine them explicitly
!-------------------------------------------------------------------------------
  use mod_scale_merge, only  :  dminfo, &
                                get_scale_xy_rd, &
                                get_scale_dm, &
                                destroy_scale_dm, &
                                r_sngl, r_dble
  implicit none

  integer :: PRC_NUM_X, PRC_NUM_Y 
  integer :: KMAX, IMAX, JMAX
  integer :: IHALO, JHALO

  integer :: nvrs2d   ! number of single precision real var
  integer :: nvrd2d   ! number of double precision real var
  !character(10),allocatable :: name_rs2d(:), name_rd2d(:)
  character(15),allocatable :: name_rs2d(:), name_rd2d(:)
  real(r_dble), allocatable :: v2dg(:,:,:)  ! nxg*nyg*(nrdv2d)
  real(r_dble), allocatable :: rd2d(:,:)    ! tmp
  real(r_sngl), allocatable :: rs2d(:,:)    ! tmp

  integer :: i, k
  integer :: luin = 10
  integer :: luout = 20
  character(10) :: ftype            ! topo/init/landuse/init_bdy/history
  character(30) :: fnin, fnout, fnconf

  namelist /dm/ PRC_NUM_X, PRC_NUM_Y, KMAX, IMAX, JMAX, nvrd2d, &
                IHALO, JHALO, fnin, fnout, ftype, &
                nvrs2d ! not used because SCALE doesn't have 2d float vars
                       ! reserved only for compatibility with old namelist
  namelist /var_real_single/ name_rs2d
  namelist /var_real_double/ name_rd2d


!-------------------------------------------------------------------------------
!-----0. configuration
  !fnin = "topo.peXXXXXX.nc"; ftype  ="landuse"; fnout ="scale_whole.grd"
  fnin  =""; ftype  =""; fnout ="scale_whole.grd"
  PRC_NUM_X = -1; PRC_NUM_Y = -1
  IHALO = -1; JHALO = -1; KMAX = -1; IMAX = -1; JMAX = -1
  nvrd2d = 0 ! topo, lon, lat
  nvrs2d = 0

  if ( command_argument_count()<1 ) then
     write(6,*) "[usage] ./merge_sub.x merge.conf"
     stop 1
  endif
  call get_command_argument(1,fnconf)
  print*, "fnconf=", trim(fnconf)

  open(luin,file=trim(fnconf),action="read")
  read(luin,nml=dm)
  if (PRC_NUM_X<=0 .or. PRC_NUM_Y<=0 .or. IHALO<0 .or. JHALO<0 .or. &
      KMAX<=0 .or. IMAX<=0 .or. JMAX <=0 .or. nvrd2d<=0 .or. &
      len_trim(ftype)==0 ) then
     write(6,*) "[error] merge: inappropriate settings in config list"
     stop 1
  else
     nvrd2d=max(nvrd2d,0)
     fnin=trim(ftype)//".peXXXXXX.nc"
     write(6,nml=dm)
  endif

  if (nvrd2d>0) then 
     allocate(name_rd2d(nvrd2d))
     rewind(luin); read(luin,nml=var_real_double)
     write(6,nml=var_real_double)
  endif

!-------------------------------------------------------------------------------
!-----1. fill the whole domains with the subdomains
  allocate( v2dg(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y,nvrd2d) )
  v2dg = 0.d0

  call get_scale_dm( dminfo, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, IHALO, JHALO )
  if (nvrd2d>0) then
      allocate( rd2d(size(v2dg,1),size(v2dg,2)) )
      do k =1, nvrd2d
         call get_scale_xy_rd( dminfo, ftype, name_rd2d(k), rd2d )
         v2dg(:,:,k) = rd2d
      enddo
      deallocate( rd2d )
  endif    

!-------------------------------------------------------------------------------
!-----2. print some diagnositic files
  do k=1, nvrd2d
     write(6,*) trim(name_rd2d(k)), ": min, max=", &
                minval(v2dg(:,:,k)), maxval(v2dg(:,:,k))
  enddo

!-------------------------------------------------------------------------------
!-----3. write into grads files
  allocate( rs2d(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y) )
  open(luout,file=trim(fnout),action="write",form="unformatted",access="direct", &
       recl=IMAX*PRC_NUM_X*JMAX*PRC_NUM_Y*r_sngl )
  i=0
  do k=1, nvrd2d
     i=i+1
     rs2d(:,:) = real( v2dg(:,:,k), r_sngl )
     write(luout,rec=i) rs2d
  enddo

  close(luout)

  call destroy_scale_dm( dminfo )
  deallocate( v2dg, rs2d )
  if (allocated(name_rs2d)) deallocate(name_rs2d)
  if (allocated(name_rd2d)) deallocate(name_rd2d)
 

endprogram
