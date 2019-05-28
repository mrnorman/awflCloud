# The AWFL Cloud Dynamical Core

A stratified, non-hydrostatic, fully compressible 3-D Cartesian grid cloud model dynamical core using [A]DER [W]ENO [F]inite-Vo[L]ume (AWFL) Numerics. The numerical scheme is designed to perform well on accelerators by performing significant amounts of local work in between DRAM memory accesses. It is written in C++ using the Kokkos library to accelerate the computations on a variety of hardware backends.
