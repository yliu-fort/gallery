Sample: fluidsGL
Minimum spec: SM 3.0

An example of fluid simulation using CUDA and CUFFT, with OpenGL rendering.

Key concepts:
Graphics Interop
CUFFT Library
Physically-Based Simulation

###
Copyright 2019 Yuxuan Liu
Modified fluidsGL_kernels.cu, implemented unconditionally stable 2nd-order advection scheme.

Selle A, Fedkiw R, Kim B, et al. An unconditionally stable MacCormack method[J]. Journal of Scientific Computing, 2008, 35(2-3): 350-371.
