#!/bin/bash 


. ./config.pps1.sh

srclist="simple_merge.f90 simple_merge_sub.f90 simple_merge_hist.f90 merge_hist_xy_2d.f90 merge_init_zxy_3d.f90"


for fn in $srclist ;
do
srcname=$fn
exename="$(echo $srcname| cut -d "." -f1).x"

cmd="$FC $FCOPTS $FCFLAGS  \
     mod_nc.f90 mod_scale_merge.f90 $srcname \
     -o $exename \
     $LDFLAGS"
echo $cmd
$cmd
echo "--------------------------------------------------------------------"

done
