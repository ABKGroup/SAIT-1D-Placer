#!/bin/bash
#source ~/SCRIPT/source_openroad.sh
export TESTCASE=1
# export MOVE=$2
# /home/memzfs_projects/SAIT_DPO/sakundu/1D_SA_GH/OpenROAD/build/src/openroad run_openroadSA.tcl -log ./log/run_test_${TESTCASE}_xxM_move.log
/home/memzfs_projects/SAIT_DPO/zhiang/SA1D/OpenROAD/build/src/openroad run_openroadSA.tcl -log ./log/run_test_${TESTCASE}_xxM_move.log
