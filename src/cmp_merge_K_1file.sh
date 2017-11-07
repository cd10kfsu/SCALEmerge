#!/bin/bash 

if [ $# -ne 1 ]; then
   echo "[usage]: $0 srcname"
   echo "if srcname is more than one files, use cmp_merge_K.sh instead"
   exit 1
fi
srcname=$1
exename="$(echo $srcname| cut -d "." -f1).x"

. ./config.pps1.sh

cmd="$FC $FCOPTS $FCFLAGS  \
     mod_nc.f90 mod_scale_merge.f90 $srcname \
     -o $exename \
     $LDFLAGS"
echo $cmd
$cmd


