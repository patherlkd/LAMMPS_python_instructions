# LAMMPS_python_instructions (see relevant pdf)
Simple instructions of how I downloaded and installed python interface for LAMMPS.

There is also an example python code for a basic LAMMPS simulation using the Python interface.

## pylammpsmpi

I have added a MWE of a simple LJ fluid working with the pylammpsmpi new interface: https://github.com/pyiron/pylammpsmpi/tree/main?tab=readme-ov-file

An example of a resulting configuration is added to the folder.

For the pylammpsmpi example you need a working mpi4py and lammps (for python) installed and working correctly. Eg. you should be able to run (in a python script or terminal):
```
from lammps import lammps
from mpi4py import MPI
```

You also need the pylammpsmpi package installed. You can so using pip:

```
pip install -U pylammpsmpi
```
The example also can be used to cupy (a GPU accelerate numpy) whereby numpy commands can be exactly replaced by cupy:

```
import cupy as cp

a = cp.array([1.,2.])
```
For example.
