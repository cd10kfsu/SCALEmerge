module mod_scale_merge
!$$$
!          .               .               .
! programer: da,cheng    org: umd    date:2017-8-17
! 
! purpose:
!   modules providing utilities to merge all scale subdomains
!   together. 
!     
!           N
!      -----------            PRC_NUM_Y-1
!      |         |
! W    |         |   E
!      |         |
!      -----------             0
!          S
!     0-------->PRC_NUM_X-1


!
!$$$
  use mod_nc, only : nc_readVar2d_real4, &
                     nc_readVar3d_real4, &
                     nc_readVar4d_real4, &
                     nc_readVar2d_real8, &
                     nc_readVar3d_real8, &
                     nc_readVar4d_real8

  implicit none


  private

  public :: r_dble, r_sngl
  public :: t_dmnode

  public :: dminfo
  public :: get_scale_dm
  public :: destroy_scale_dm
  public :: get_scale_xyt_rs
  public :: get_scale_xy_rd
  public :: get_scale_zxy_rd


  integer,parameter :: r_dble = kind(0.d0)
  integer,parameter :: r_sngl = kind(0.e0)
  integer,parameter :: r_size = r_dble
  integer,parameter :: nvmax = 20

  type :: t_dmnode
    logical :: linit = .false.
    integer :: np    = 0
    integer :: imax  = 0
    integer :: jmax  = 0
    integer :: kmax  = 0
    integer :: prc_num_x = 0
    integer :: prc_num_y = 0
    logical,allocatable :: prc_s(:), prc_n(:), prc_w(:), prc_e(:)
    integer,allocatable :: nxl(:), nyl(:)
    integer,allocatable :: isx(:), iex(:), jsx(:), jex(:)  
    integer,allocatable :: isg(:), ieg(:), jsg(:), jeg(:)
    integer,allocatable :: PRC_2Drank(:,:)
  endtype

  type(t_dmnode),save :: dminfo

  interface destroy_scale_dm
    module procedure destroy_dmnode
  endinterface


contains 

!-------------------------------------------------------------------------------
! double SFLX_SW_dif(y, x)
!-------------------------------------------------------------------------------
subroutine get_scale_xy_rd( d, ftype, varname, v2d )
 implicit none

  type(t_dmnode),intent(in)    :: d
  character(*),  intent(in)    :: ftype
  character(*),  intent(in)    :: varname
  real(r_dble),  intent(inout) :: v2d(:,:)  !xy

  character(200) :: fnin
  integer :: ipes, ipee 
  integer :: i
  integer :: isg, ieg, jsg, jeg, isx, iex, jsx, jex
  real(r_dble),allocatable :: tmpv2d(:,:)

  if ( size(v2d,1) /= d%prc_num_x*d%imax .or. &
       size(v2d,2) /= d%prc_num_y*d%jmax ) then
     write(*,*) "[error] get_scale_xy_rd: dim of input var /= size of SCALE domain"
     stop 1
  endif

  v2d = 0.d0
  fnin=trim(ftype)//".peXXXXXX.nc"
  ipes=len_trim(ftype)+4; ipee=ipes+5
  do i = 0, dminfo%np - 1

     write(fnin(ipes:ipee),"(I6.6)") i

     isg = dminfo%isg(i); ieg = dminfo%ieg(i); jsg = dminfo%jsg(i); jeg = dminfo%jeg(i)
     isx = dminfo%isx(i); iex = dminfo%iex(i); jsx = dminfo%jsx(i); jex = dminfo%jex(i)

     allocate( tmpv2d(d%nxl(i),d%nyl(i)) )
     tmpv2d = 0.d0
     call nc_readVar2d_real8(  trim(fnin), trim(varname), tmpv2d )
     v2d(isg:ieg,jsg:jeg) = tmpv2d(isx:iex,jsx:jex)
     deallocate( tmpv2d )

  enddo

endsubroutine

!-------------------------------------------------------------------------------
! read vars such as float QS(time, z, y, x)
!-------------------------------------------------------------------------------
subroutine get_scale_xyzt_rs( d, ftype, varname, v4d )
 implicit none

  type(t_dmnode),intent(in)    :: d
  character(*),  intent(in)    :: ftype
  character(*),  intent(in)    :: varname
  real(r_sngl),  intent(inout) :: v4d(:,:,:,:)  !xyzt

  character(200) :: fnin
  integer :: ipes, ipee 
  integer :: i
  integer :: isg, ieg, jsg, jeg, isx, iex, jsx, jex
  real(r_sngl),allocatable :: tmpv4d(:,:,:,:)

  if ( size(v4d,1) /= d%prc_num_x*d%imax .or. &
       size(v4d,2) /= d%prc_num_y*d%jmax .or. &
       size(v4d,3) /= d%kmax ) then
     write(*,*) "[error] get_scale_xyt_rs: dim of input var /= size of SCALE domain"
     stop 1
  endif

  v4d = 0.d0
  fnin=trim(ftype)//".peXXXXXX.nc"
  ipes=len_trim(ftype)+4; ipee=ipes+5
  do i = 0, dminfo%np - 1

     write(fnin(ipes:ipee),"(I6.6)") i
     isg = dminfo%isg(i); ieg = dminfo%ieg(i); jsg = dminfo%jsg(i); jeg = dminfo%jeg(i)
     isx = dminfo%isx(i); iex = dminfo%iex(i); jsx = dminfo%jsx(i); jex = dminfo%jex(i)
     allocate( tmpv4d(d%nxl(i),d%nyl(i),size(v4d,3),size(v4d,4)) )
     tmpv4d = 0.d0
     call nc_readVar4d_real4(  trim(fnin), trim(varname), tmpv4d )
     v4d(isg:ieg,jsg:jeg,:,:) = tmpv4d(isx:iex,jsx:jex,:,:)
     deallocate( tmpv4d )

  enddo

endsubroutine


!-------------------------------------------------------------------------------
! read vars like MSLP(time, y, x)
!-------------------------------------------------------------------------------
subroutine get_scale_xyt_rs( d, ftype, varname, v3d )
 implicit none

  type(t_dmnode),intent(in)    :: d
  character(*),  intent(in)    :: ftype
  character(*),  intent(in)    :: varname
  real(r_sngl),  intent(inout) :: v3d(:,:,:)  !xyt

  character(200) :: fnin
  integer :: ipes, ipee 
  integer :: i
  integer :: isg, ieg, jsg, jeg, isx, iex, jsx, jex
  real(r_sngl),allocatable :: tmpv3d(:,:,:)  !xyt

  if ( size(v3d,1) /= d%prc_num_x*d%imax .or. &
       size(v3d,2) /= d%prc_num_y*d%jmax ) then
     write(*,*) "[error] get_scale_xyt_rs: dim of input var /= size of SCALE domain"
     stop 1
  endif

  v3d = 0.d0
  fnin=trim(ftype)//".peXXXXXX.nc"
  ipes=len_trim(ftype)+4; ipee=ipes+5
  do i = 0, dminfo%np - 1

     write(fnin(ipes:ipee),"(I6.6)") i

     isg = dminfo%isg(i); ieg = dminfo%ieg(i); jsg = dminfo%jsg(i); jeg = dminfo%jeg(i)
     isx = dminfo%isx(i); iex = dminfo%iex(i); jsx = dminfo%jsx(i); jex = dminfo%jex(i)

     allocate( tmpv3d(d%nxl(i),d%nyl(i),size(v3d,3)) )
     tmpv3d = 0.d0
     call nc_readVar3d_real4(  trim(fnin), trim(varname), tmpv3d )
     v3d(isg:ieg,jsg:jeg,:) = tmpv3d(isx:iex,jsx:jex,:)
     deallocate( tmpv3d )

  enddo

endsubroutine


!-------------------------------------------------------------------------------
! read var like double DENS(y, x, z)
! will output vars as (x,y,z)
!-------------------------------------------------------------------------------
subroutine get_scale_zxy_rd( d, ftype, varname, v3d )
 implicit none

  type(t_dmnode),intent(in)    :: d
  character(*),  intent(in)    :: ftype
  character(*),  intent(in)    :: varname
  real(r_dble),  intent(inout) :: v3d(:,:,:)  ! xyz

  character(200) :: fnin
  integer :: ipes, ipee 
  integer :: i, k 
  integer :: isg, ieg, jsg, jeg, isx, iex, jsx, jex
  real(r_dble),allocatable :: tmpv3d(:,:,:)   ! zxy

  if ( size(v3d,1) /= d%prc_num_x*d%imax .or. &
       size(v3d,2) /= d%prc_num_y*d%jmax ) then
     write(*,*) "[error] get_scale_xyt_rs: dim of input var /= size of SCALE domain"
     stop 1
  endif

  v3d = 0.d0
  fnin=trim(ftype)//".peXXXXXX.nc"
  ipes=len_trim(ftype)+4; ipee=ipes+5
  do i = 0, dminfo%np - 1

     write(fnin(ipes:ipee),"(I6.6)") i

     isg = dminfo%isg(i); ieg = dminfo%ieg(i); jsg = dminfo%jsg(i); jeg = dminfo%jeg(i)
     isx = dminfo%isx(i); iex = dminfo%iex(i); jsx = dminfo%jsx(i); jex = dminfo%jex(i)

     allocate( tmpv3d(size(v3d,3),d%nxl(i),d%nyl(i)) )
     tmpv3d = 0.d0
     call nc_readVar3d_real8(  trim(fnin), trim(varname), tmpv3d )
     do k = 1, size(tmpv3d,1)
        v3d(isg:ieg,jsg:jeg,k) = tmpv3d(k,isx:iex,jsx:jex)
     enddo
     deallocate( tmpv3d )

  enddo

endsubroutine




!-------------------------------------------------------------------------------
! return the handler <d> carring the grid mapping infos from the whole domain
! to the subdomains
!-------------------------------------------------------------------------------
subroutine get_scale_dm( d, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, IHALO, JHALO, KMAX )
 implicit none

  type(t_dmnode),intent(inout) :: d
  integer,intent(in) :: PRC_NUM_X, PRC_NUM_Y
  integer,intent(in) :: IMAX, JMAX
  integer,intent(in) :: IHALO, JHALO
  integer,intent(in),optional :: KMAX

  character(30) :: fnin  ="topo.peXXXXXX.nc"

  integer ::  i

  d%prc_num_x = PRC_NUM_X
  d%prc_num_y = PRC_NUM_Y
  d%imax      = IMAX
  d%jmax      = JMAX
  if (present(KMAX)) d%kmax = KMAX
  d%np = PRC_NUM_X*PRC_NUM_Y
  !print*, "nproc =", d%np
  call init_dmnode( d )

!-----1.1 calculate 2d rank for each process
! copied from the scale_rm_process.F90 /scalelib/src/atmos-rm/communication
  do i = 0, d%np-1

     d%PRC_2Drank(i,1) = mod(i,PRC_NUM_X)
     d%PRC_2Drank(i,2) = (i-d%PRC_2Drank(i,1))/PRC_NUM_X
     write(6,*) "i, prc_2drank=", i, d%PRC_2Drank(i,1), d%PRC_2Drank(i,2)

  enddo

     
!           N
!      -----------            PRC_NUM_Y-1
!      |         |
! W    |         |   E
!      |         |
!      -----------             0
!          S
!     0-------->PRC_NUM_X-1

  do i = 0, d%np-1

     write(fnin(8:13),"(I6.6)") i
     !write(6,*) "read file: ", TRIM(fnin)

     d%isx(i)=1; d%iex(i) = IMAX; d%jsx(i)=1; d%jex(i) = JMAX 
     !-----1.2 determine the size of subdomain ( include halo or not ) from the
     !         netcdf file.
     !         also the 
     if ( d%PRC_2Drank(i,1) == 0 .or. d%PRC_2Drank(i,1) == PRC_NUM_X-1 ) then
        d%nxl(i) = IMAX + IHALO
        if ( d%PRC_2Drank(i,1) == 0 ) then
           d%prc_w(i) = .true.
           d%isx(i) = 1+IHALO; d%iex(i) = IMAX+IHALO
        else
           d%prc_e(i) = .true.
           d%isx(i) = 1      ; d%iex(i) = IMAX
        endif
     else
        d%nxl(i) = IMAX
     endif

     if ( d%PRC_2Drank(i,2) == 0 .or. d%PRC_2Drank(i,2) == PRC_NUM_Y-1 ) then
        d%nyl(i) = JMAX + JHALO
        if ( d%PRC_2Drank(i,2) == 0 ) then
           d%prc_s(i) = .true.
           d%jsx(i) = 1+JHALO; d%jex(i) = JMAX+JHALO
        else
           d%prc_n(i) = .true.
           d%jsx(i) = 1      ; d%jex(i) = JMAX
        endif
     else
        d%nyl(i) = JMAX
     endif

     !write(6,"(A,A,I2,A,I2)") TRIM(fnin), ": x=", nxl, ", y=", nyl

     !-----1.3 determine the locations of the subdomains in the whole domain
     d%isg(i) = 1    + d%PRC_2Drank(i,1)*IMAX
     d%ieg(i) = IMAX + d%PRC_2Drank(i,1)*IMAX
     d%jsg(i) = 1    + d%PRC_2Drank(i,2)*JMAX
     d%jeg(i) = JMAX + d%PRC_2Drank(i,2)*JMAX
     write(6,*) "-----------------------------------------------------------------"
     !write(6,*) "proc_1d, proc_2d_i, proc_2d_j = ", i, PRC_2Drank(i,1), PRC_2Drank(i,2)
     write(6,*) "proc_1d = ", i
     write(6,"(2A,4(I4))") "fnin, isg, ieg, jsg, jeg = ", TRIM(fnin), d%isg(i), d%ieg(i), d%jsg(i), d%jeg(i)
     write(6,"(2A,4(L4))") "fnin,   w,   e,   s,   n = ", TRIM(fnin), d%prc_w(i), d%prc_e(i), d%prc_s(i), d%prc_n(i)
     write(6,"(2A,4(I4))") "fnin, isx, iex, jsx, jex = ", TRIM(fnin), d%isx(i), d%iex(i), d%jsx(i), d%jex(i)

  enddo

endsubroutine

 

subroutine init_dmnode( d )
  implicit none

  type(t_dmnode),intent(inout) :: d

  if ( d%linit ) then
     write(6,*) "[warning] init_dnmode: dm_node already initialized. return..."
     return
  endif

  allocate( d%PRC_2Drank(-1:d%np-1,2), &
            d%prc_w(0:d%np-1), &
            d%prc_e(0:d%np-1), &
            d%prc_s(0:d%np-1), &
            d%prc_n(0:d%np-1), &
            d%isx(0:d%np-1), &
            d%iex(0:d%np-1), &
            d%jsx(0:d%np-1), &
            d%jex(0:d%np-1), &
            d%isg(0:d%np-1), & 
            d%ieg(0:d%np-1), &
            d%jsg(0:d%np-1), &
            d%jeg(0:d%np-1), &
            d%nxl(0:d%np-1), &
            d%nyl(0:d%np-1)    )

  d%PRC_2Drank(:,:) = -1
  d%prc_w = .false. ; d%prc_e = .false. ; d%prc_s = .false. ; d%prc_n = .false. 
  d%isx = -1; d%iex = -1; d%jsx = -1; d%jex = -1
  d%isg = -1; d%ieg = -1; d%jsg = -1; d%jeg = -1
 
  d%linit = .true.

endsubroutine


subroutine destroy_dmnode( d )
   implicit none

   type(t_dmnode),intent(inout) :: d

  if ( .not. d%linit ) then
     write(6,*) "[warning] init_dnmode: dm_node not initialized yet. return...."
     return
  endif

  deallocate( d%PRC_2Drank, &
              d%prc_w, d%prc_e, d%prc_s, d%prc_n, &
              d%isx, d%iex, d%jsx, d%jex, &
              d%isg, d%ieg, d%jsg, d%jeg, &
              d%nxl, d%nyl &
            )
 d%linit = .false.

endsubroutine


endmodule
  

 

