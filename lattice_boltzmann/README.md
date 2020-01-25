Copyright 2019 Yuxuan Liu  
Implemented lattice boltzmann method with modern opengl.  
Dimension and lattice set implemented: D2Q9,D3Q27  
Boundary conditions implemented: periodic, pressure drop periodic, noslip, inlet, reflective outlet, charcteristic outlet.  
Stability enhancement implemented: explicit entropy correction (support high Reynolds flow simulation)  
Curved boundary: staircase  
Corner: haven't been treated specifically, detect&increase viscosity in 3d to stabilize simulation in extremely high Re  
Boundary method: link-wise  
No-slip boundary method: bounce-back (2nd-order accurary has been validated in matlab code)  

## matlab codes  
run system("nvcc -ptx <cu_name>.cu") to compile ptx files before running the script.  
make sure that nvidia compiler has been installed.  
