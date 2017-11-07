program main
! dirty code: wrapper to merge all SCALE subdomain together
  use mod_nc,          only : nc_readVar2d_real8, &
                              nc_readVar2d_real4
  use mod_scale_merge 
  implicit none

  integer :: PRC_NUM_X, PRC_NUM_Y 
  integer :: KMAX, IMAX, JMAX
  integer :: IHALO, JHALO

  integer :: nvrs2d   ! number of single precision real var
  integer :: nvrd2d   ! number of double precision real var
  !character(10),allocatable :: name_rs2d(:), name_rd2d(:)
  character(15),allocatable :: name_rs2d(:), name_rd2d(:)
  real(r_dble), allocatable :: v2dg(:,:,:)  ! nxg*nyg*(nrdv2d+nvrs2d)
  real(r_dble), allocatable :: rd2d(:,:)    ! tmp
  real(r_sngl), allocatable :: rs2d(:,:)    ! tmp

  integer :: isx, iex, jsx, jex, isg, ieg, jsg, jeg
  integer :: ipes, ipee
  integer :: i, k
  integer :: luin = 10
  integer :: luout = 20
  character(10) :: ftype            ! topo/init/landuse/init_bdy/history
  character(30) :: fnin, fnout, fnconf


  namelist /dm/ PRC_NUM_X, PRC_NUM_Y, KMAX, IMAX, JMAX, nvrs2d, nvrd2d, &
                IHALO, JHALO, fnin, fnout, ftype
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
     write(6,*) "[usage] ./merge.x merge.conf"
     stop 1
  endif
  call get_command_argument(1,fnconf)
  print*, "fnconf=", trim(fnconf)

  open(luin,file=trim(fnconf),action="read")
  read(luin,nml=dm)
  if (PRC_NUM_X<=0 .or. PRC_NUM_Y<=0 .or. IHALO<0 .or. JHALO<0 .or. &
      KMAX<=0 .or. IMAX<=0 .or. JMAX <=0 .or. (nvrs2d+nvrd2d)<=0 .or. &
      len_trim(ftype)==0 ) then
     write(6,*) "[error] merge: inappropriate settings in config list"
     stop 1
  else
     nvrd2d=max(nvrd2d,0)
     nvrs2d=max(nvrs2d,0)
     fnin=trim(ftype)//".peXXXXXX.nc"
     write(6,nml=dm)
  endif

  if (nvrd2d>0) then 
     allocate(name_rd2d(nvrd2d))
     rewind(luin); read(luin,nml=var_real_double)
     write(6,nml=var_real_double)
  endif

  if (nvrs2d>0) then
     allocate(name_rs2d(nvrs2d))
     rewind(luin); read(luin,nml=var_real_single)
     write(6,nml=var_real_single)
  endif
  close(luin)

!-------------------------------------------------------------------------------
!-----1. fill the whole domains with the subdomains
  allocate( v2dg(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y,nvrs2d+nvrd2d) )
  v2dg = 0.d0

  call get_scale_dm( dminfo, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, IHALO, JHALO )

  ipes=len_trim(ftype)+4; ipee=ipes+5
  do i = 0, dminfo%np - 1

     write(fnin(ipes:ipee),"(I6.6)") i

     isg = dminfo%isg(i); ieg = dminfo%ieg(i); jsg = dminfo%jsg(i); jeg = dminfo%jeg(i)
     isx = dminfo%isx(i); iex = dminfo%iex(i); jsx = dminfo%jsx(i); jex = dminfo%jex(i)

     if (nvrd2d>0) then
        allocate( rd2d(dminfo%nxl(i),dminfo%nyl(i)) )
        rd2d = 0.d0
        do k = 1, nvrd2d
           call nc_readVar2d_real8( trim(fnin), trim(name_rd2d(k)), rd2d )
           v2dg(isg:ieg,jsg:jeg,k) = rd2d(isx:iex,jsx:jex)
        enddo
        deallocate( rd2d )
     endif

     if (nvrs2d>0) then
        allocate( rs2d(dminfo%nxl(i),dminfo%nyl(i)) )
        rs2d = 0.d0
        do k = 1, nvrs2d
           call nc_readVar2d_real4( trim(fnin), trim(name_rs2d(k)), rs2d )
           v2dg(isg:ieg,jsg:jeg,nvrd2d+k) = rs2d(isx:iex,jsx:jex)
        enddo
        deallocate( rs2d )
     endif

  enddo
 
!-------------------------------------------------------------------------------
!-----2. print some diagnositic files
  do k=1, nvrd2d
     write(6,*) trim(name_rd2d(k)), ": min, max=", &
                minval(v2dg(:,:,k)), maxval(v2dg(:,:,k))
  enddo

  do k=1, nvrs2d
     write(6,*) trim(name_rs2d(k)), ": min, max=", &
                minval(v2dg(:,:,nvrd2d+k)), maxval(v2dg(:,:,nvrd2d+k))
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
  do k=1, nvrs2d
     i=i+1
     rs2d(:,:) = real( v2dg(:,:,nvrd2d+k), r_sngl )
     write(luout,rec=i) rs2d
  enddo

  close(luout)

  call destroy_scale_dm( dminfo )
  deallocate( v2dg, rs2d )
  if (allocated(name_rs2d)) deallocate(name_rs2d)
  if (allocated(name_rd2d)) deallocate(name_rd2d)
 

endprogram
