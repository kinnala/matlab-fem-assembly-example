# Triangular FEM assembly in 2D
Example of 2D Poisson equation in Matlab/Octave.

I have used this for teaching and student projects, some ideas I've used with students:
- Assemble linear elasticity 2x2 block system
- Assemble Stokes 3x3 block system using the unstable P1-P1 element and find its nullspace
- Assemble Stokes 3x3 block system using P1-P1 stab. (Brezzi-Pitkäranta element) and solve lid driven cavity
- Extend code to 3D assembly and solve 3D Poisson problem
- Extend code to quadratic triangular element, e.g., allow switching P2 for u and/or v to use Taylor-Hood for Stokes block system

Playing with short codes like this eventually led me to create scikit-fem which supports all sorts of elements and bilinear forms: https://github.com/kinnala/scikit-fem
