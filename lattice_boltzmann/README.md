Copyright 2019 Yuxuan Liu  
Implemented lattice boltzmann method with modern opengl.  
Dimension and lattice set implemented: D2Q9,D3Q27  
Boundary conditions implemented: periodic, pressure drop periodic, noslip, inlet, reflective outlet, charcteristic outlet.  
Stability enhancement implemented: explicit entropy correction (support high Reynolds flow simulation)  
Curved boundary: staircase  
Corner: haven't been treated specifically, detect&increase viscosity in 3d to stabilize simulation in extremely high Re  
Boundary method: link-wise  
No-slip boundary method: bounce-back (2nd-order accurary has been validated in matlab code)  

## Compile
First make sure you have CMake, Git, and GCC by typing as root (sudo) apt-get install g++ cmake git  
and then get the required packages: Using root (sudo) and type apt-get install libsoil-dev libglm-dev libassimp-dev libglew-dev libglfw3-dev libxinerama-dev libxcursor-dev libxi-dev  
