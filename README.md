This codes use fewest switches surface hopping method ( Molecular dynamics with electronic friction, John C. Tully, 1990) to simulate polariton dynamics where electronic excitation of two site is
coupled directly with the photon mode. 

This code is written in fortran90. Python is used only for parallel running. For running the code, use makefile. Then run aout file. 
Necessary changes can be made depending the system in the input file (AFSSH.inp) and in the hamiltonian of the file mod_afssh.f90

