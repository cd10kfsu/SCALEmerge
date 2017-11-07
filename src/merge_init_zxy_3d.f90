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

  real(r_dble), allocatable :: v3dg(:,:,:,:) ! nxg*nyg*nz*nv
  real(r_dble), allocatable :: tmpv3d(:,:,:) ! nxg*nyg*nz
  real(r_sngl), allocatable :: buf2d(:,:) ! nxg*nyg

  integer :: i, k, n
  integer :: luout = 20
  character(10) :: ftype            ! topo/init/landuse/init_bdy/boundary/history
  character(30) :: fnout

  namelist /dm/ PRC_NUM_X, PRC_NUM_Y, KMAX, IMAX, JMAX, &
                IHALO, JHALO, fnout, ftype, varname

!-------------------------------------------------------------------------------
!-----0. configuration
  ftype     = "init"
  fnout     = "scale_whole_Q.grd"
  PRC_NUM_X = 4
  PRC_NUM_Y = 3
  IHALO     = 2
  JHALO     = 2 
  KMAX      = 36
  IMAX      = 60
  JMAX      = 60

  nv        = 6
  varname(1)= "DENS"
  varname(2)= "QC"
  varname(3)= "QR"
  varname(4)= "QI"
  varname(5)= "QS"
  varname(6)= "QG"
 
  write(*,nml=dm)
 
!-------------------------------------------------------------------------------
!-----1. fill the whole domains with the subdomains
  allocate(    v3dg(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y,KMAX,nv) )
  allocate(  tmpv3d(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y,KMAX)    )
  v3dg = 0.d0

  call get_scale_dm( dminfo, PRC_NUM_X, PRC_NUM_Y, IMAX, JMAX, IHALO, JHALO, KMAX )
  do n = 1, nv
     call get_scale_zxy_rd( dminfo, ftype, varname(n), tmpv3d )
     v3dg(:,:,:,n) = tmpv3d
  enddo

!-------------------------------------------------------------------------------
!-----3. write into grads files
  allocate( buf2d(IMAX*PRC_NUM_X,JMAX*PRC_NUM_Y) )
  open(luout,file=trim(fnout),action="write",form="unformatted",access="direct", &
       recl=IMAX*PRC_NUM_X*JMAX*PRC_NUM_Y*r_sngl )
  print*, "kmax=", dminfo%kmax
  i=0
  do n=1, nv
     do k = 1, dminfo%kmax
        i=i+1
        buf2d(:,:) = real( v3dg(:,:,k,n), r_sngl )
        write(luout,rec=i) buf2d
     enddo
  enddo

  close(luout)

endprogram
