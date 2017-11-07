program main
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
  use mod_scale_merge 
  implicit none

  integer :: PRC_NUM_X, PRC_NUM_Y 
  integer :: KMAX, IMAX, JMAX
  integer :: IHALO, JHALO

  integer :: nv
  character(15) :: varname(10)

  real(r_dble), allocatable :: v2dg(:,:,:) ! nxg*nyg*nv
  real(r_dble), allocatable :: tmpv2d(:,:) ! nxg*nyg
  real(r_sngl), allocatable :: buf2d(:,:) ! nxg*nyg

  integer :: i, n
  integer :: luout = 20
  character(10) :: ftype            ! topo/init/landuse/init_bdy/boundary/history
  character(30) :: fnout

  namelist /dm/ PRC_NUM_X, PRC_NUM_Y, KMAX, IMAX, JMAX, &
                IHALO, JHALO, fnout, ftype, nv, varname

!-------------------------------------------------------------------------------
!-----0. configuration
  ftype     = "history"
  fnout     = "scale_whole.grd"
  PRC_NUM_X = 4
  PRC_NUM_Y = 3
  IHALO     = 0
  JHALO     = 0 
  KMAX      = 36
  IMAX      = 60
  JMAX      = 60

  nv        = 4
  varname(1)= "topo"
  varname(2)= "lon"
  varname(3)= "lat"
  varname(4)= "lsmask"

  write(*,nml=dm)
 
!-------------------------------------------------------------------------------
!-----1. fill the whole domains with the subdomains
  allocate(    v2dg(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y,nv) )
  allocate(  tmpv2d(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y)    )
  v2dg = 0.d0

  call get_scale_dm( dminfo, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, IHALO, JHALO )
  do n = 1, nv
     call get_scale_xy_rd( dminfo, ftype, varname(n), tmpv2d )
     v2dg(:,:,n) = tmpv2d
  enddo

!-------------------------------------------------------------------------------
!-----3. write into grads files
  allocate( buf2d(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y) )
  open(luout,file=trim(fnout),action="write",form="unformatted",access="direct", &
       recl=IMAX*PRC_NUM_X*JMAX*PRC_NUM_Y*r_sngl )
  i=0
  do n=1, nv
     i=i+1
     buf2d(:,:) = real( v2dg(:,:,n), r_sngl )
     write(luout,rec=i) buf2d
  enddo

  close(luout)

endprogram
