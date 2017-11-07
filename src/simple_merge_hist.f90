program main
! dirty code: wrapper to merge all SCALE subdomain together
  use mod_nc,          only : nc_readVar3d_real8, &
                              nc_readVar3d_real4
  use mod_scale_merge 
  implicit none

  integer :: PRC_NUM_X, PRC_NUM_Y 
  integer :: KMAX, IMAX, JMAX
  integer :: IHALO, JHALO
  integer :: nt

  integer :: nvrs2dd   ! number of single precision real var
  integer :: nvrd2dd   ! number of double precision real var
  !character(10),allocatable :: name_rs2d(:), name_rd2d(:)
  character(15),allocatable :: name_rs2dd(:), name_rd2dd(:)
  real(r_dble), allocatable :: v2ddg(:,:,:,:)  ! nxg*nyg*nt*(nrdv2d+nvrs2d)
  real(r_dble), allocatable :: rd2dd(:,:,:)    ! tmp
  real(r_sngl), allocatable :: rs2dd(:,:,:)    ! tmp
  real(r_sngl), allocatable :: rsbuf2d(:,:)

  integer :: isx, iex, jsx, jex, isg, ieg, jsg, jeg
  integer :: ipes, ipee
  integer :: i, k, n
  integer :: luin = 10
  integer :: luout = 20
  character(10) :: ftype            ! topo/init/landuse/init_bdy/history
  character(30) :: fnin, fnout, fnconf


  namelist /dm/ PRC_NUM_X, PRC_NUM_Y, KMAX, IMAX, JMAX, nvrs2dd, nvrd2dd, &
                IHALO, JHALO, fnin, fnout, ftype, nt
  namelist /hist_var_real_single/ name_rs2dd
  namelist /hist_var_real_double/ name_rd2dd


!-------------------------------------------------------------------------------
!-----0. configuration
  !fnin = "topo.peXXXXXX.nc"; ftype  ="landuse"; fnout ="scale_whole.grd"
  fnin  =""; ftype  =""; fnout ="scale_whole.grd"
  PRC_NUM_X = -1; PRC_NUM_Y = -1
  IHALO = -1; JHALO = -1; KMAX = -1; IMAX = -1; JMAX = -1
  nvrd2dd = 0 ! topo, lon, lat
  nvrs2dd = 0

  if ( command_argument_count()<1 ) then
     write(6,*) "[usage] ./merge.x merge.conf"
     stop 1
  endif
  call get_command_argument(1,fnconf)
  print*, "fnconf=", trim(fnconf)

  open(luin,file=trim(fnconf),action="read")
  read(luin,nml=dm)
  if (PRC_NUM_X<=0 .or. PRC_NUM_Y<=0 .or. IHALO<0 .or. JHALO<0 .or. &
      KMAX<=0 .or. IMAX<=0 .or. JMAX <=0 .or. (nvrs2dd+nvrd2dd)<=0 .or. &
      nt==0 .or. len_trim(ftype)==0 ) then
     write(6,*) "[error] merge: inappropriate settings in config list"
     stop 1
  else
     nvrd2dd=max(nvrd2dd,0)
     nvrs2dd=max(nvrs2dd,0)
     fnin=trim(ftype)//".peXXXXXX.nc"
     write(6,nml=dm)
  endif

  if (nvrd2dd>0) then 
     allocate(name_rd2dd(nvrd2dd))
     rewind(luin); read(luin,nml=hist_var_real_double)
     write(6,nml=hist_var_real_double)
  endif

  if (nvrs2dd>0) then
     allocate(name_rs2dd(nvrs2dd))
     rewind(luin); read(luin,nml=hist_var_real_single)
     write(6,nml=hist_var_real_single)
  endif
  close(luin)

!-------------------------------------------------------------------------------
!-----1. fill the whole domains with the subdomains
  allocate( v2ddg(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y,nt,nvrs2dd+nvrd2dd) )
  v2ddg = 0.d0

  !call get_scale_dm( dminfo, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, IHALO, JHALO )
  call get_scale_dm( dminfo, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, 0, 0 )

  ipes=len_trim(ftype)+4; ipee=ipes+5
  do i = 0, dminfo%np - 1

     write(fnin(ipes:ipee),"(I6.6)") i

     isg = dminfo%isg(i); ieg = dminfo%ieg(i); jsg = dminfo%jsg(i); jeg = dminfo%jeg(i)
     isx = dminfo%isx(i); iex = dminfo%iex(i); jsx = dminfo%jsx(i); jex = dminfo%jex(i)

     if (nvrd2dd>0) then
        allocate( rd2dd(dminfo%nxl(i),dminfo%nyl(i),nt) )
        rd2dd = 0.d0
        do k = 1, nvrd2dd
           call nc_readVar3d_real8( trim(fnin), trim(name_rd2dd(k)), rd2dd )
           v2ddg(isg:ieg,jsg:jeg,:,k) = rd2dd(isx:iex,jsx:jex,:)
        enddo
        deallocate( rd2dd )
     endif

     if (nvrs2dd>0) then
        allocate( rs2dd(dminfo%nxl(i),dminfo%nyl(i),nt) )
        rs2dd = 0.d0
        do k = 1, nvrs2dd
           call nc_readVar3d_real4( trim(fnin), trim(name_rs2dd(k)), rs2dd )
           v2ddg(isg:ieg,jsg:jeg,:,nvrd2dd+k) = rs2dd(isx:iex,jsx:jex,:)
        enddo
        deallocate( rs2dd )
     endif

  enddo
 
!-------------------------------------------------------------------------------
!-----2. print some diagnositic files
  do n=1, nt
     do k=1, nvrd2dd
        write(6,*) trim(name_rd2dd(k)), ": min, max=", &
                minval(v2ddg(:,:,n,k)), maxval(v2ddg(:,:,n,k))
     enddo
  enddo
  do n=1, nt
     do k=1, nvrs2dd
        write(6,*) trim(name_rs2dd(k)), ": min, max=", &
                minval(v2ddg(:,:,n,nvrd2dd+k)), maxval(v2ddg(:,:,n,nvrd2dd+k))
     enddo
  enddo

!-------------------------------------------------------------------------------
!-----3. write into grads files
  allocate( rsbuf2d(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y) )
  open(luout,file=trim(fnout),action="write",form="unformatted",access="direct", &
       recl=IMAX*PRC_NUM_X*JMAX*PRC_NUM_Y*r_sngl )
  i=0
  do n=1, nt
     do k=1, nvrd2dd+nvrs2dd
        i=i+1
        rsbuf2d(:,:) = real( v2ddg(:,:,n,k), r_sngl )
        write(luout,rec=i) rsbuf2d
     !enddo
     !do k=1, nvrs2dd
     !   i=i+1
     !   rsbuf2d(:,:) = real( v2ddg(:,:,n,nvrd2dd+k), r_sngl )
     !   write(luout,rec=i) rsbuf2d
     enddo
  enddo
  close(luout)

  call destroy_scale_dm( dminfo )
  deallocate( v2ddg, rsbuf2d )
  if (allocated(name_rs2dd)) deallocate(name_rs2dd)
  if (allocated(name_rd2dd)) deallocate(name_rd2dd)
 

endprogram
