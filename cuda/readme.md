## Several practices for gpu/cuda programming.

### 1. Spectral-projection method
Original code from cuda-10.0 code example.   
Original code applied method of charcteristics + projection, results in a super viscous fluid.  
Modifications have been made to implement an unconditionally stable, 2nd-order advection scheme.  

### 2. nbody
Original code from cuda-10.0 code example.   
Original code applied brute-force algorithm to compute two-body interactions.  
Modifications have been made to implement an hierachical method to approximate n-body gravity and reduce time complexity from O(n^2) to O(nlogn).  
