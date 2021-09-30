#!/bin/bash -e

#--- compile unit test code
make gtest

declare EXE_DIR='./test_bin'
export OMP_NUM_THREADS=2

#--- PS::Vector3<T> extension
${EXE_DIR}/gtest_vec_ext

#--- COMM_TOOL::
${EXE_DIR}/gtest_comm_tool_SerDes
mpirun -np 2 ${EXE_DIR}/gtest_comm_tool_broadcast
mpirun -np 2 ${EXE_DIR}/gtest_comm_tool_allGather

#--- MD_EXT::boltzmann_dist
${EXE_DIR}/gtest_blz_dist

#--- MD_EXT::CellIndex
${EXE_DIR}/gtest_cell_index

#--- static array for FullParticle
${EXE_DIR}/gtest_fixed_vector
${EXE_DIR}/gtest_basic_connect

#--- IntraPair::
mpirun -np 2 ${EXE_DIR}/gtest_intra_pair

#--- force test
mpirun -np 2 ${EXE_DIR}/gtest_force_LJ
mpirun -np 2 ${EXE_DIR}/gtest_force_coulomb
mpirun -np 2 ${EXE_DIR}/gtest_force_bond
mpirun -np 2 ${EXE_DIR}/gtest_force_angle
mpirun -np 2 ${EXE_DIR}/gtest_force_dihedral
mpirun -np 2 ${EXE_DIR}/gtest_force_improper

#--- file I/O test
mpirun -np 2 ${EXE_DIR}/gtest_fileIO
