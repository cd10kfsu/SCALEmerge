#!/bin/bash 

nodename=`uname -n`
#--------------------------------------------------
# K pps1
#--------------------------------------------------

NCDIR="/volume64/data/ra001011/cda/work/pkg/netcdf-4.4.1.1-intel/"
H5DIR="/volume64/data/ra001011/cda/work/pkg/hdf5-1.8.19-intel"
ZDIR="/volume64/data/ra001011/cda/work/pkg/zlib-1.2.11-intel"

NC_FCFLAG="-I$NCDIR/include"
H5_FCFLAG="-I$H5DIR/include"
Z_FCFLAG="-I$ZDIR/include"

NC_LDFLAG="-L$NCDIR/lib -lnetcdff -lnetcdf"
H5_LDFLAG="-L$H5DIR/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5"
Z_LDFLAG="-L$ZDIR/lib -lz"

FCFLAGS="$Z_FCFLAG  \
         $H5_FCFLAG \
         $NC_FCFLAG"

LDFLAGS="$NC_LDFLAG \
         $H5_LDFLAG \
         $Z_LDFLAG"

FC="ifort"
FCOPTS="-assume byterecl -check bounds -warn all"


echo "--------------------------------------------------------------------"
echo "            config lists for node: $nodename"
echo 
echo "FCFLAGS=$FCFLAGS"
echo 
echo "LDFLAGS=$LDFLAGS"
echo 
echo "FC=$FC"
echo 
echo "FCOPTS=$FCOPTS"
echo 
echo "--------------------------------------------------------------------"

# one example below

#cmd="$FC $FCOPTS $FCFLAGS  mod_nc.f90 mod_scale_merge.f90 simple_merge_hist.f90 -o merge_hist_K.x $LDFLAGS"
#echo $cmd
#$cmd


