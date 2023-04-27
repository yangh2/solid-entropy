#!/bin/bash


cwd=`pwd`

rm */{Eigenvalues.dat,covariance_pair.dat,density_matrix.dat,density_matrix_row1.dat,force_matrix.dat,xdat-ent.out,Trun,Mass}

cd bcc_Li
echo "Running examples in bcc_Li"
../../bin/xdat-ent-init.sh > /tmp/NUL
../../bin/xdat-ent-rhoALL_quantum_mod -f xyz.ideal XDATCAR -s symmetry_operation.dat > xdat-ent.out
cd $cwd

cd fcc_Al
echo "Running examples in fcc_Al"
../../bin/xdat-ent-init.sh > /tmp/NUL
../../bin/xdat-ent-rhoALL_quantum_mod -f xyz.ideal XDATCAR -s symmetry_operation.dat > xdat-ent.out
sleep 3;
../../bin/freq2quan.py > freq.out;
../../bin/force2dos xyz.ideal 24 > dos.out;
../../bin/force2band xyz.ideal > band.out
cd $cwd

cd hcp_Ti
echo "Running examples in hcp_Ti"
../../bin/xdat-ent-init.sh > /tmp/NUL
../../bin/xdat-ent-rhoALL_quantum_mod -f xyz.ideal XDATCAR -s symmetry_operation.dat > xdat-ent.out
cd $cwd

cd HEA_MoNbTaW
echo "Running examples in HEA_MoNbTaW"
../../bin/xdat-ent-init.sh > /tmp/NUL
sleep 1;
../../bin/xdat-ent-rhoALL_quantum_alloy -f xyz.ideal XDATCAR -s symmetry_operation.dat  > xdat-ent.out
cd $cwd
