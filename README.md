# LAMMPS_python_instructions
Simple instructions of how I downloaded and installed python interface for LAMMPS.

### RELEVANT FOR MYRIAD, change as needed

0. Ensure you have wget

1. Download lammps locally

wget https://github.com/lammps/lammps/archive/refs/heads/develop.zip

2. Send lammps to your HPC (or skip if using on your local machine)

scp develop.zip <hpc-address>:/home/<username>

3. login to your HPC and unzip the file in a directory called lammps/

cd lammps/
mkdir lammps-develop-python
unzip develop.zip -d ./lammps-develop-python/

4. cd into the src/ folder in lammps-develop-python

cd lammps-develop-python/lammps-develop/src/

5. include python library (exact path for Python.h) in the c compiler for lammps. open src/MAKE/Makefile.mpi. Change the CCFLAGS line like so: 

CCFLAGS =       -g -O3 -std=c++11 -I/lustre/shared/ucl/apps/python/3.11.4/gnu-4.9.2/include/python3.11 -lpython3.11

6. make lammps in shared mode with the python package (still within src/ folder)

machine="mpi" 
echo " == Making LAMMPS-$machine for $nodename == "                                                                                                                                                                
echo "Machine filename: $machine"                                                                                                                                                                                  
make clean-${machine} 
make no-all                                                                                                                                                                                                        
make yes-EXTRA-FIX                                                                                                                                                                                                 
make yes-EXTRA-PAIR                                                                                                                                                                                                
make yes-MANIFOLD                                                                                                                                                                                                  
make yes-MOLECULE                                                                                                                                                                                                  
make yes-RIGID                                                                                                                                                                                                     
make yes-MANYBODY                                                                                                                                                                                                  
make yes-BROWNIAN                                                                                                                                                                                                  
make yes-PYTHON
make -j 4 mode=shared ${machine}                                                   
                                                                                                                                             
make install-python

7. Run a test lammps python from one of the examples located in lammps-develop/python/examples or any examples given to you.
