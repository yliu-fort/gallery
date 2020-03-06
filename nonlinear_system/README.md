Copyright 2020 Yuxuan Liu  
Implemented 4/5th Adaptive Runge-kutta Integration  
Integration performed on GPU may diverge heavily, be cautious  
that a 16-loop upper limit is applied in kernel function  
in order to suppress the branch divergence thus there is no  
guarantee for accuracy/stability if large timestep or tolerance is applied  
uncomment "#define VISUAL" in main.cpp in order to turn on real-time rendering  
