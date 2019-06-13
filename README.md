# The AWFL Cloud Dynamical Core

## Note that this is a work in progress

A stratified, non-hydrostatic, fully compressible 3-D Cartesian grid cloud model dynamical core using [A]DER [W]ENO [F]inite-Vo[L]ume (AWFL) Numerics. The numerical scheme is designed to perform well on accelerators by performing significant amounts of local work in between DRAM memory accesses. It is written in C++ using the Kokkos library to accelerate the computations on a variety of hardware backends.

Some attributes of the numerical scheme:
* __ADER__: ADER time stepping using Differential Transforms for single-step, single-stage high-order time stepping
* __WENO__: A new function-based Weighted Essentially Non-Oscillatory (WENO) limiting for crisp resolution of discontinuities without oscillations. 
* __Finite-Volume__: Flux-vector Finite-Volume spatial integration for local conservation
* __Upwind__ fluxes using the upwind Godunov linear state via flux-vector splitting for slight dissipation in the numerical scheme
* __High-Order__: Up to ninth-order accuracy with WENO limiting and 15th-order accuracy without WENO limiting
* __Compressible__: Solves the fully compressible Euler equations for applicability to the full range of atmospheric motions
* __Non-Hydrostatic__: Resolves non-hydrostatic motions, thus admitting vertically propagating acoustic waves
* __Dimensionally-Split__ for reduced cost
* Uses __Perturbed__ quantities with the hydrostatic background state removed for better resolution of hydrostatic balance
* Includes a __2-D__ option to simulate the traditional NH atmospheric test cases
* Uses __MPI__ (nearest neighbor communication) for off-node parallelism
* Uses __Kokkos__ for on-node parallelism
