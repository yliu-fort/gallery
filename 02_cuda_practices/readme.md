## Several practices for gpu/cuda programming.

### 1. Spectral-projection method
Original code from cuda-10.0 code example.   
Original code applied method of charcteristics + projection, results to super viscous behaviour.  
Modifications have been made to implement an unconditionally stable, 2nd-order advection scheme.  

### 2. nbody
Original code from cuda-10.0 code example.   
Original code applied brute-force algorithm to compute two-body interactions.  
Modifications have been made to implement an octree method with parallel construction&traversal to approximate n-body gravity and reduce time complexity from O(n^2) to O(nlogn).  
