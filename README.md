## SCALEmerge

#### name: Da,Cheng (cda@umd.edu)    org: umd

A simple utility to combine all SCALE subdomain fields to 
the whole-domain fields purely in Fortran. No MPI needed.


Tested with SCALE TOPO LANDUSE, HISTORY, INIT files. Only
works for the vars with coordinates related to x, y, z

#### 0. Structure:
```
src      # source codes
vars     # vars in various SCALE files
testcase # a test case to check if the program runs correctly
```

#### 1. Compile:
  ```
cd src
vim config.pps1.sh      # modify config.pps1.sh
./cmp_merge_K_all.sh    # compiles all the utilities
  ```
  
#### 2. Run the test case:
  ```
cd testcase
./simple_merge_sub.x merge.conf_topo
#then use grads to open valid.ctl & new.ctl, check if vars in new.grd is the same as the those in valid.grd
  ```


