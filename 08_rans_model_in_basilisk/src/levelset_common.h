#ifndef LEVELSET_COMMON_H
#define LEVELSET_COMMON_H

#define BICUBIC 1
#include "alimare/alex_functions.h"
#include "alimare/LS_reinit.h"
// distance to interface
// d = alpha*Delta/sqrt(sq(n.x) + sq(n.y) + sq(n.z))

#define INTERFACE_CELL_DISTANCE (1.0)
#define INTERFACE_DAMPING_FACTOR (10.0)
#define INTERFACE_CELL_SHIFT (1.0)

struct LevelSetConstants {
	double icl; // interface cell refinement level, cell width = s/(1<<icl)
	double icw;
};

struct LevelSetConstants ls_constant = {
.icl = 0, // interface cell width, 1/Delta, for interface turbulence damping
.icw = 1
};

scalar dist[];
face vector dist_f[];

#ifdef TWO_PHASE
void compute_distance_to_interface() {
  if (ls_constant.icl == 0) {
    foreach()
      dist[] = 1.0;
     foreach_face()
      dist_f.x[] = 1.0;
  } else {
    foreach(){
      dist[] = f[] - 0.5;
    }
    //boundary({dist});

    int nbit = LS_reinit(dist,it_max=20);

    foreach_face() {
      dist_f.x[] = fabs(0.5 * (dist[] + dist[-1]));
      double d = fabs(dist_f.x[]) / ls_constant.icw;
      dist_f.x[] = d < INTERFACE_CELL_SHIFT ? 1.0 : 
      ((exp(-(INTERFACE_DAMPING_FACTOR/INTERFACE_CELL_DISTANCE)*sq(d-INTERFACE_CELL_SHIFT)) - exp(-INTERFACE_DAMPING_FACTOR))/(1.0 - exp(-INTERFACE_DAMPING_FACTOR)));
      dist_f.x[] = max(dist_f.x[], 0.0);
      dist_f.x[] = min(dist_f.x[], 1.0);
    }

    foreach() {
      double d = fabs(dist[]) / ls_constant.icw;
      dist[] = d < INTERFACE_CELL_SHIFT ? 1.0 : 
      ((exp(-(INTERFACE_DAMPING_FACTOR/INTERFACE_CELL_DISTANCE)*sq(d-INTERFACE_CELL_SHIFT)) - exp(-INTERFACE_DAMPING_FACTOR))/(1.0 - exp(-INTERFACE_DAMPING_FACTOR)));
      dist[] = max(dist[], 0.0);
      dist[] = min(dist[], 1.0);
    }
  }
    //boundary({dist});
    //boundary({dist_f});
}
#else
void compute_distance_to_interface() {
  foreach()
    dist[] = 0;
  foreach_face()
    foreach_dimension()
      dist_f.x[] = 0;
}
#endif

event init (i=0,last) {
	// k, omega should be initialized by user
    compute_distance_to_interface();
}

#if TREE
// mesh adaptation happens before solving viscous term
event adapt (i++) {
  compute_distance_to_interface();
}
#else
event compute_distance (i++, last) {
  compute_distance_to_interface();
}
#endif

#endif