# The AWFL Cloud Dynamical Core

## Note that this is a work in progress

A stratified, non-hydrostatic, fully compressible 3-D Cartesian grid cloud model dynamical core using [A]DER [W]ENO [F]inite-Vo[L]ume (AWFL) Numerics. The numerical scheme is designed to perform well on accelerators by performing significant amounts of local work in between DRAM memory accesses. It is written in C++ using the [YAKL](github.com/mrnorman/YAKL) portability library to accelerate the computations on Nvidia and AMD GPU backends.

Some attributes of the numerical scheme:
* __ADER__: ADER time stepping using Differential Transforms for single-step, single-stage high-order time stepping with a maximum stable CFL value of one. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2017MS001247
* __WENO__: A new function-based Weighted Essentially Non-Oscillatory (WENO) limiting for crisp resolution of discontinuities without oscillations. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2017MS001247
* __Finite-Volume__: Flux-vector Finite-Volume spatial integration for local conservation to machine precision
* __Upwind__ fluxes using the upwind Godunov linear state via a characteristics-based flux-vector splitting that gives enough dissipation to render stability for any spatial order of accuracy even on unstaggered grids.
* __High-Order__: Up to ninth-order spatial accuracy with WENO limiting
* __Compressible__: Solves the fully compressible Euler equations for applicability to the full range of atmospheric motions
* __Non-Hydrostatic__: Resolves non-hydrostatic motions, thus admitting vertically propagating acoustic waves
* __Dimensionally-Split__ for reduced cost with minimal impact on accuracy due to Cartesian spatial coordinates
* Uses __Perturbed__ quantities with the hydrostatic background state removed for better resolution of hydrostatic balance
* Includes a __2-D__ option to simulate the traditional NH atmospheric test cases
* Uses __MPI__ (nearest neighbor communication) for off-node parallelism
* Uses __C++ Portability__ for on-node parallelism

## Future Work
The most pressing developments needed are:
1) A more efficient handling of fast acoustics, particularly in the vertical direction. This could include implicit, multi-level, and Lagrangian treatments.
2) A variable vertical grid spacing (simply a matter of altering the SageMath scripts to compute Vandermonde matrices that respond to grid spacing and inverting those at model initialization).
3) Add tracer transport routines
4) Create a plan for incorporating moisture, possibly through an equivalent potential temperature.

## Software Dependencies
* Nvidia GPUs require CUB: https://github.com/NVlabs/cub
* AMD GPUs require hipCub and rocPrim:
  * https://github.com/ROCmSoftwarePlatform/hipcub
  * https://github.com/ROCmSoftwarePlatform/rocPRIM

## How to Build
```
ln -s build/[whatever_machine].inc mach.inc
make
```

## How to build on Summit
```bash
git clone https://github.com/mrnorman/awflCloud.git
cd awfulCloud
git submodule update --init
module load parallel-netcdf gcc cuda
ln -s build/summit_gpu.inc mach.inc
## edit const.h and change "ord" to 3,5,7, or 9
## edit input.txt and change doWeno to true or false
make clean
make -j
jsrun -n 1 -a 1 -c 1 -g 1 ./cloudFV
```

* Requires an MPI library (that provides mpic++) as well as parallel-netcdf (https://trac.mcs.anl.gov/projects/parallel-netcdf)

