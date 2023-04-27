#!/bin/bash

cd src;
make clean;
make
cp xdat-ent-rhoALL_quantum_mod xdat-ent-rhoALL_quantum_npt xdat-ent-rhoALL_quantum_alloy aflow-sym.sh force2dos force2band freq2quan.py xdat-ent-init.sh ../bin/
cd -;
