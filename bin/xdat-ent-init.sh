#!/bin/bash

chtab="
#ATOMIC_NUM MASS DIAM (A) SHORT_NAME LONG_NAME ALIAS\n
  0    1          0       Vac Vacancy\n
  1    1.007900   0.741    H Hydrogen ?\n
  2    4.002600   3.000   He Helium ?\n
  3    6.941000   3.039   Li Lithium l\n
  4    9.012180   2.226   Be Beryllium S\n
  5   10.810000   1.799    B Boron s   # ORIGINAL: 1.589\n
  6   12.011000   1.426    C Carbon s  \n
  7   14.006740   1.098    N Nitrogen s\n
  8   15.999400   1.5    O Oxygen s 1.207    O Oxygen s\n
  9   18.998403   1.418    F Fluorine s\n
 10   20.179000   3.000   Ne Neon ?\n
 11   22.989770   3.716   Na Sodium L\n
 12   24.305000   3.197   Mg Magnesium l\n
 13   26.981540   2.863   Al Aluminum M\n
 14   28.086000   2.352   Si Silicon S \n
 15   30.973760   2.210    P Phosphorus S\n
 16   32.060000   2.050    S Sulfur s\n
 17   35.453000   1.891   Cl Chlorine s\n
 18   39.948000   3.000   Ar Argon ?\n
 19   39.098000   4.544    K Potassium ?\n
 20   40.080000   3.947   Ca Calcium L\n
 21   44.955900   3.212   Sc Scandium l\n
 22   47.900000   2.896   Ti Titanium M\n
 23   50.941400   2.622    V Vanadium m\n
 24   51.996000   2.498   Cr Chromium m\n
 25   54.938000   2.731   Mn Manganese M\n
 26   55.847000   2.482   Fe Iron m\n
 27   58.933200   2.506   Co Cobalt m\n
 28   58.700000   2.492   Ni Nickel m\n
 29   63.546000   2.556   Cu Copper m\n
 30   65.380000   2.665   Zn Zinc m\n
 31   69.720000   3.025   Ga Gallium M modified for atomic volume\n
 32   72.590000   2.450   Ge Germanium\n
 33   74.921600   2.490   As Arsenic ?\n
 34   78.960000   2.321   Se Selenium ?\n
 35   79.904000   2.284   Br Bromine ?\n
 36   83.800000   4.040   Kr Krypton ?\n
 37   85.467800   4.700   Rb Rubidium ?\n
 38   87.620000   4.100   Sr Strontium ?\n
 39   88.905900   3.551    Y Yttrium L\n
 40   91.220000   3.179   Zr Zirconium l\n
 41   92.906400   2.858   Nb Niobium M\n
 42   95.940000   2.725   Mo Molybdenum M\n
 43   97.000000   2.703   Tc Technetium ?\n
 44  101.070000   2.650   Ru Ruthenium M\n
 45  102.905500   2.690   Rh Rhodium M\n
 46  106.400000   2.751   Pd Palladium M\n
 47  107.868000   2.889   Ag Silver M\n
 48  112.400000   2.979   Cd Cadmium M\n
 49  114.820000   3.251   In Indium l\n
 50  118.690000   2.810   Sn Tin M\n
 51  121.750000   2.900   Sb Antimony M\n
 52  127.600000   2.860   Te Tellurium ?\n
 53  126.904500   3.000    I Iodine ?\n
 54  131.300000   3.000   Xe Xenon ?\n
 55  132.905400   5.309   Cs Cesium H\n
 56  137.340000   4.347   Ba Barium H\n
 57  138.905500   3.739   La Lanthanum ?\n
 58  140.120000   3.650   Ce Cerium ?\n
 59  140.907700   3.640   Pr Praseodynium ?\n
 60  144.240000   3.628   Nd Neodymium L\n
 61  145.000000   3.000   Pm Promethium ?\n
 62  150.400000   3.579   Sm Samarium L\n
 63  151.960000   3.989   Eu Europium L\n
 64  157.250000   3.573   Gd Gadolinium L\n
 65  158.925400   3.525   Tb Terbium L\n
 66  162.500000   3.503   Dy Dysprosium L\n
 67  164.930400   3.486   Ho Holmium L\n
 68  167.260000   3.468   Er Erbium L\n
 69  168.934200   3.447   Tm Thulium L\n
 70  173.040000   3.880   Yb Ytterbium L\n
 71  174.970000   3.435   Lu Lutetium ?\n
 72  178.490000   3.127   Hf Hafnium ?\n
 73  180.947900   2.860   Ta Tantalum ?\n
 74  183.500000   2.741    W Tungsten M\n
 75  186.207000   2.741   Re Rhenium M\n
 76  190.200000   2.675   Os Osmium M\n
 77  192.220000   2.714   Ir Iridium M\n
 78  195.090000   2.775   Pt Platinum M\n
 79  196.966500   2.884   Au Gold M\n
 80  200.590000   3.005   Hg Mercury l\n
 81  204.370000   3.408   Tl Thallium L\n
 82  207.200000   3.500   Pb Lead L\n
 83  208.980400   3.090   Bi Bismuth ?\n
 84  209.000000   3.345   Po Polonium ?\n
 85  210.000000   3.000   At Astatine ?\n
 86  222.000000   3.000   Rn Radon ?\n
 87  223.000000   3.000   Fr Francium ?\n
 88  226.025400   4.458   Ra Radium ?\n
 89  227.000000   3.756   Ac Actinium\n
 90  232.038100   3.595   Th Thorium ?\n
 91  231.035900   3.212   Pa Protactinium ?\n
 92  238.029000   2.770    U Uranium ?\n
 93  237.048200   2.620   Np Neptunium ?\n
 94  244.000000   3.026   Pu Plutonium ?\n
 95  243.000000   3.450   Am Americium ?\n
 96  247.000000   3.477   Cm Curium ?\n
 97  247.000000   3.398   Bk Berkelium ?\n
 98  251.000000   3.377   Cf Californium ?\n
 99  252.000000   3.000   Es Einsteinium ?\n
100  257.000000   3.000   Fm _Fermium ?\n
101  258.000000   3.000   Md _Mendelevium ?\n
102  259.000000   3.000   No _Nobelium ?\n
103  262.000000   3.000   Lr _Lawrencium ?\n
104  261.000000   3.000   Rf _Rutherfordium ?\n
105  262.000000   3.000   Db _Dubnium ?\n
106  263.000000   3.000   XX _Seaborgium ?\n
107  262.000000   3.000   Bf _Bohrium ?\n
108  265.000000   3.000   Hs _Hassium ?\n
109  265.000000   3.000   Mt _Meitnerium ?\n
110    1.000001          3.0     Xm Mix ?\n
"

if [ -f Trun ];
then
    echo "Trun file already exists!"
fi

if ! [ -f Trun ];
then
    echo "Trun file has been generated!"
    grep TEBEG OUTCAR | awk '{print $3}' | sed 's/;//' > Trun;
fi

if [ -f Mass ];
then
    echo "Mass file already exists!"
fi
#echo -e $chtab
if ! [ -f Mass ];
then
    echo "Mass file has been generated!"
    ntype=`awk 'NR==6{print $0}' POSCAR | wc -w`;
    echo $ntype > Mass;
    for i in `awk 'NR==6{print $0}' POSCAR`;
    do
	#echo $i;
	an=`echo -e $chtab | awk -v ch=$i '$4==ch{print $1}'`;
	echo -n "$an " >> Mass
    done
    echo >> Mass;
    grep POMASS OUTCAR | awk 'BEGIN{FS="= "}; {print $2}' | sed 's/....../& /g' >> Mass
fi

echo "Before you run:"
echo "  Check lattice constants and ideal sites in xyzfile"
echo '  Check symmetry-file (directional space). See "aflow-sym.sh"'
echo ;
echo "========================================"
echo "Run with:";
echo "  xdat-ent-rhoALL_quantum_mod -f xyzfile -s symmetry-file XDATCAR* > out"
echo "========================================"
echo ;
echo "Be careful with:"
echo "  Use xdat-ent-rhoALL_quantum_alloy, if chemical species swapping is allowed"
echo "  Use xdat-ent-rhoALL_quantum_npt for npt simulation"
echo ;
echo 'Also see "freq2quan.py","E0EK.sh" '


