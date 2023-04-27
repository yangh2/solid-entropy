           __________________________________________________

             VIBRATIONAL ENTROPY OF CRYSTALLINE SOLIDS FROM
                   COVARIANCE OF ATOMIC DISPLACEMENTS

                               Yang Huang
           __________________________________________________


                            <2023-01-18 Wed>


Table of Contents
_________________

1. Purpose
2. Contents
3. Required packages
4. How to install
5. List of tools
6. Pay attention to





1 Purpose
=========

  This archive contains tools to calculate entropy from covariance
  matrix based on method describe in Huang, Y.; Widom, M. Vibrational
  Entropy of Crystalline Solids from Covariance of Atomic Displacements.
  Entropy 2022, 24, 618. <https://doi.org/10.3390/e24050618>. These
  tools are designed especially for output from MD simulation
  implemented in VASP. You are free to modify and redistribute the code
  for general purpose.


2 Contents
==========

  Following list of files should be contained in the directory:
  - README

  - install.sh

  - src

    Contain source files.

  - example

    Contain examples

  - bin

    Contain compiled binaries.


3 Required packages
===================

  1) GSL - GNU Scientific Library

     The GNU Scientific Library (GSL) is a numerical library for C and
     C++ programmers. It is free software under the GNU General Public
     License. It is installed on most linux systems. To install it,
     check <https://www.gnu.org/software/gsl/>.

  2) Aflow [optional]

     The aflow-sym.sh requires AFLOW codes in order to run properly. To
     install it, check <http://aflow.org/install-aflow/>.

  3) pos2xyz and xyz2pos from euler

     The xdat-ent program is originally written on euler so that it
     utilizes certain mgtools on euler. Unfortunately the author is
     unwilling to modify these codes to make it more convenient for
     general users and you will have to ask for permissions to access
     these tools.

     The "pos2xyz" converts a POSCAR file to a XYZ file. The "xyz2pos"
     converts a XYZ file to a POSCAR file.


4 How to install
================

  Go to the main directory and run "./install.sh". It will compile the
  source files and put binaries into the bin directory.

  After you run the command, you run the test script "./test.sh" in the
  examples directory.


5 List of tools
===============

  [BEFORE YOU RUN]
  - xdat-ent-init.sh

    - Description

      Preparing "Trun" and "Mass" input files for
      "xdat-ent-rhoALL_quantum_mod" commands.

    - Input:
      - INCAR
      - POTCAR
      - /home/Tools/src/mgtools/INFO/CHTAB.

    - Usage:

      xdat-ent-init.sh *Only work on euler.

    - Output:

      Generate "Trun" and "Mass".

  - aflow-sym.sh
    - Description

      This script utilizes "aflow" tools to help find symmetry
      operations. It is recommended to run on a small cell for
      efficiency.

    - Input:

      POSCAR or XYZ file.

    - Usage:

      aflow-sym.sh [ -x <xyz file> ] [ -p <POSCAR file> ]

    - Output:

      A "symmetry_operation.dat " which contains number of operations N
      and N 3x3 operation matrix in direction space.

  [WHEN YOU RUN]
  - xdat-ent-rhoALL_quantum_mod tools

    - Description

      + xdat-ent-rhoALL_quantum_mod

        Pure element crystalline from NVT simulation

      + xdat-ent-rhoALL_quantum_npt

        Pure element crystalline from NPT simulation

      + xdat-ent-rhoALL_quantum_alloy

        Alloys from NVT simulation.

      xdat-ent-rhoALL_quantum_mod has been tested to be robust while
      only a few examples have been performed with
      xdat-ent-rhoALL_quantum_npt and xdat-ent-rhoALL_quantum_alloy
      command. So be careful with these two commands!

    - Input:

      "example" directory provides examples of input files.

      - Trun

        This file contains in value which is the temperature of MD
        simulation.

      - Mass

        This file contains three rows:
        1. type of species.
        2. atomic numbers.
        3. atomic masses [a.u.] in atomic units.

      - xyzfile

        This file contains the ideal structure of a crystal. The first
        three rows are lattice constants. The 4th row contains number of
        atoms N. The following N rows should obeys: dx,dy,dz,na ...
        where dx,dy,dz is directional vector and na is atomic number.

      - symmetry-file

        This file contains symmetry operation matrices U. The first row
        is number of operations N and then it is followed by N 3x3
        matrices U each of which is a operation matrix in directional
        space (b'=Ub).

      - XDATCAR

        Output from VASP MD simulations.


  - Usage:

    xdat-ent-rhoALL_quantum_mod -f xyzfile -s symmetry-file XDATCAR* >
    out

  - Output:

    + covariance_pair.dat

      3x3 covariance matrices in real space of all symmetry
      distinguishing pairs.

    + density_matrix.dat

      See reference for more mathematical description. A 3Nx3N matrix.

    + force_matrix.dat

      See reference for more mathematical description. A 3Nx3N force
      matrix.

    + Eigenvalues.dat

      A 3Nx2 matrix. The first row are eigenvalues of density matrix.
      The second row are frequencies (hbar*omega [eV]) of eigen modes
      See reference for mathematical relation.

    + "out":

      The last line of "out" reports calculated entropy with S_classic
      (classical entropy) and S_quantum (quantum entropy)

  [AFTER YOU RUN]
  - freq2quan.py

    - Description

      Calculate free energy F[eV/at], internal energy U [eV/at], entropy
      S [kB/at] and heat capacity Cv [kB/at] from Eigenvalues.dat

    - Input: "Trun", "Eigenvalues.dat"

    - Usage: freq2quan.py

    - Output: calculated quantities.

  - force2dos

    - Description

      Calculate free energy F[eV/at], internal energy U [eV/at], entropy
      S [kB/at] and heat capacity Cv [kB/at] from Eigenvalues.dat

    - Input:

      - "Trun", "force_matrix.dat", xyzfile

      - KDim (an integer indicate kmesh KDimxKDimxKdim).

      - recip_vectors of a primitive cell

        3x3 matrix

    - Usage:

      force2dos xyzfile KDim > out;

    - Output:

      In "out" file, there are two columns of datas. The first column is
      hbar*omega [meV]. The second one is vibrational DOS
      [states/atom/meV].

  - force2band

    - Description

      Calculate free energy F[eV/at], internal energy U [eV/at], entropy
      S [kB/at] and heat capacity Cv [kB/at] from Eigenvalues.dat

    - Input:

      - "Trun", "force_matrix.dat", xyzfile

      - KDim (an integer indicate kmesh KDimxKDimxKdim).

      - recip_vectors of a primitive cell

        3x3 matrix

      - kpath

    - Usage:

      force2band xyzfile > out;

    - Output:

      Eigenvalues along kpath.


6 Pay attention to
==================

  - symmetry operations:

    Symmetry operations will greatly improve precision of calculated
    entropy for a short run at a cost of calculation time. With symmetry
    operation turned on, each configuration will be process N times
    where N is the number of operations.

    This program is known to give incorrect results when applied to
    low-symmetry super cell from a highly symmetrical crystalline
    structure, e.g. a 5x5x5 primitive cell of BCC Li. Hence, it is
    recommended to run simulation at its high-symmetry unit cell. e.g. a
    super cell of a 2-atom BCC unit cell or a 4-atom FCC unit cell.

  - Eigenvalues.dat:

    A reasonable calculation should give 3 zero eigenvalues in the first
    column in Eigenvalues.dat and correspondingly 3 nan or large number
    in the second column. So be careful with zero or nan values in
    Eigenvalues.dat.

    Imaginary modes may show up which are denoted as negative
    eigenvalues. As long as there are 3 zero eigenvalues, we think the
    program is doing its job.
