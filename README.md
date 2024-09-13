# Meson sekelton code

This repository contains a [Meson](https://mesonbuild.com/) skeleton. It is used in the [High-Performance Computiong:  MD
with C++
project](https://pastewka.github.io/MolecularDynamics/_project/general_remarks.html).

# Compile
cd builddir
meson compile

### Run executables of parallel code, change of domain decompositions needs to be done inside the main files. Executables for all parallel code, have to be started in the folder containing xyz files, special error messages for missing files are not shown, although meson should take care of them.
cd builddir/millestones/07/
mpirun -n 4 ./milestone7
### data for plots with timestep for parallel computation on gold clusters
cd builddir/millestones/08/
mpirun -n 4 ./milestone08
### for straining, has to be started in the folder containing xyz files, special error messages for missing files are not shown
cd builddir/millestones/09/
mpirun -n 20 ./milestone09

### Plots of the resulting files are generated in results/results1.ipynb
