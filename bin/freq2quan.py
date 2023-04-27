#!/usr/bin/env python3

################################################################################
# Calculate free energy F[eV/at], internal energy U [eV/at], entropy S [kB/at] #
# and heat capacity Cv [kB/at] from Eigenvalues.dat                            #
# Usage: freq2quan.py
# Requirment: "Trun", "Eigenvalues.dat"
################################################################################

import numpy as np

kb=8.617e-5;
Trun = np.loadtxt('Trun');
beta = 1/(Trun*kb);
freq=np.loadtxt('Eigenvalues.dat');
freqev=freq[:,1];

l,= np.shape(freqev);

f = 0;
u = 0;
s = 0;
cv = 0;

for i in range(l):
    omega = freqev[i];
    if (omega > 100):
        continue;
    if (np.isnan(omega)):
        continue;
    if (omega < 0 ):
        print("Imaginary mode:", omega);
        continue;
    if (omega < 1e-2):
        continue;
    f = f+0.5*omega + 1/beta*np.log(1-np.exp(-beta*omega));
    #print (f, omega);
    u = u+0.5*omega + omega/(np.exp(beta*omega)-1);
    s = s-np.log(1-np.exp(-beta*omega)) + beta*omega/(np.exp(beta*omega)-1);
    tt=omega*beta/2
    cv = cv+(tt/np.sinh(tt))**2;
    
f = f/l*3;                      # free energy per atom [eV/at]
u = u/l*3;                      # enthalpy per atom   [eV/at]
s = s/l*3;                      # entropy per atom  [kb/at]
cv = cv/l*3;                    # heat capacity per atom [kb/at]

print("T[k] F[eV/at] U[eV/at] S[kB/at] Cv[kB/at]");
print(Trun,f,u,s,cv)
