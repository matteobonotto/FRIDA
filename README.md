# FRIDA
routines/binaries of the FRIDA code.

The FRIDA (FRee-boundary Integro-Differential Axisimmetric) code expoits a coupled Finite Element Method – Boundary Element Method (FEM-BEM) approach for the solution of the free-boundary axi-symmetric plasma equilibrium problem (i.e., the Grad-Shafranov equation).

## Meshing
FRIDA supports triagular meshes, which can be generated in several ways. However, we provide some routines to generate triangular meshes suitable for FRIDA usage. Such routines are based on the [GMSH](https://gmsh.info/) software. Up to this verison of FRIDA (v3.3), generating meshes directly via Matlab script is supported for Windows OS only. 

## Use compiled C/C++/Fortran subroutines
Some routines can be compiled into mex functions (see [here](https://it.mathworks.com/help/matlab/call-mex-functions.html?lang=en) for further details). To mex these routines, just run the matlab script `compile/compile_mex_functions.m`.

## Examples
Some examples of how to setup a geometry and run FRIDA can be found in the sub-directory `examples`. The following examples are reported:
1. static equilibrium  ...
2. time domain equilibrium ...
3. time domain vacuum run ...
 
### Citations
Please cite the following papers:
```
@article{bonotto2022coupled,
  title={A coupled FEM-BEM approach for the solution of the free-boundary axi-symmetric plasma equilibrium problem},
  author={Bonotto, M and Abate, D and Bettini, P and Villone, F and others},
  journal={Commun. Comput. Phys.},
  volume={31},
  number={1},
  pages={27--59},
  year={2022}
}

@article{bonotto2023efficient,
  title={Efficient Numerical Solution of Coupled Axisymmetric Plasma Equilibrium and Eddy Current Problems},
  author={Bonotto, Matteo and Abate, Domenico and Bettini, Paolo and Iaiunese, Antonio and Isernia, Nicola and Villone, Fabio},
  journal={IEEE Access},
  volume={11},
  pages={27489--27505},
  year={2023},
  publisher={IEEE}
}

