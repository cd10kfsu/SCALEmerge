#!/bin/bash 

. ./config.pps1.sh

cmd="$FC $FCOPTS $FCFLAGS  \
     mod_nc.f90 mod_scale_merge.f90 simple_merge.f90 \
     -o merge_K.x \
     $LDFLAGS"
echo $cmd
$cmd


