/** 
# Adapt wavelet leave interface
This is a slightly altered version of the the adapt_wavelet function used to adapt tree structures. The point of this function is to retain maximum level refinement around the interface, described by a volume fraction field vol_frac.
the padding parameter can have a integer value of 0, 1 or 2. If padding=0, maximum level will be retained only at the interface cells. If padding=1 maximum level will be retained in the interface cells and all closest neighboring cells to the interface cells. if padding=2, it will retain max level 2 cell distances away from the interface cells. You get the idea....

Modifiers: Arne Bockmann and Oystein Lande 2018
*/

struct Adapt_leave_interface {
  scalar * slist; // list of scalars
  scalar * vol_frac; // the volume fraction scalar
  double * max;   // tolerance for each scalar
  int maxlevel_air;   // relaxation factor based on vol_frac for each scalar
  int maxlevel;   // maximum level of refinement
  int minlevel;   // minimum level of refinement (default 1)
  int padding;    // number of neighbor cells to padd on each side of the interface being preserved
  scalar * list;  // list of fields to update (default all)
};


astats adapt_wavelet_leave_interface (struct Adapt_leave_interface p)
{
  if (p.list == NULL)
    p.list = all;
  if (is_constant(cm))
    restriction (p.slist);
  else {
    scalar * listr = list_concat ({cm}, p.slist);
    restriction (listr);
    free (listr);
  }


  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (!is_constant(s) && s.restriction != no_restriction)
      listc = list_add (listc, s);

  // refinement
  if (p.minlevel < 1)
    p.minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
	if (cell.flags & too_coarse) {
	  cell.flags &= ~too_coarse;
	  refine_cell (point, listc, refined, &tree->refined);
	  st.nf++;
	}
	continue;
      }
      else { // !is_leaf (cell)
	if (cell.flags & refined) {
	  // cell has already been refined, skip its children
	  cell.flags &= ~too_coarse;
	  continue;
	}
	// check whether the cell or any of its children is local
	bool local = is_local(cell);
	if (!local)
	  foreach_child()
	    if (is_local(cell)){
			local = true;
		   	break;
		}
	if (local) {
	  int i = 0;
	  static const int just_fine = 1 << (user + 3);
	  for (scalar s in p.slist) {
	    double max = p.max[i++], sc[1 << dimension];
	    int c = 0;
	    foreach_child()
	      sc[c++] = s[];
	    s.prolongation (point, s);
	    c = 0;
	    foreach_child() {
			double vf0 = 0;
			for (scalar vf in p.vol_frac) {
				vf0 = vf[];
			}
              int maxlevel = vf0 < 0.001 ? p.maxlevel_air : p.maxlevel;
	      double e = fabs(sc[c] - s[]);
	      if (e > max && level < maxlevel) {
		cell.flags &= ~too_fine;
		cell.flags |= too_coarse;
	      }
	      else if ((e <= max/1.5 || level > maxlevel) &&
		       !(cell.flags & (too_coarse|just_fine))) {
		if (level >= p.minlevel)
		  cell.flags |= too_fine;
	      }
	      else if (!(cell.flags & too_coarse)) {
		cell.flags &= ~too_fine;
		cell.flags |= just_fine;
	      }
	      // arnbo: always set interface cells to the finest level
	      for (scalar vf in p.vol_frac) {
                if (vf[] > 0.0001 && vf[] < 0.9999 && level < p.maxlevel) {
                  cell.flags |= too_coarse;
                  cell.flags &= ~too_fine;
            	  cell.flags &= ~just_fine;
                    if (p.padding > 0){
                       foreach_neighbor(p.padding){
                          cell.flags |= too_coarse;
                          cell.flags &= ~too_fine;
                          cell.flags &= ~just_fine;
                        }
                    }
                }
              }


	      s[] = sc[c++];
	    }
	  }
	  foreach_child() {
	    cell.flags &= ~just_fine;
	    if (!is_leaf(cell)) {
	      cell.flags &= ~too_coarse;
	      if (level >= p.maxlevel)
		cell.flags |= too_fine;
	    }
	    else if (!is_active(cell))
	      cell.flags &= ~too_coarse;
	  }
	}
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);
  
  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= 0; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
	if (level == l) {
	  if (!is_leaf(cell)) {
	    if (cell.flags & refined)
	      // cell was refined previously, unset the flag
	      cell.flags &= ~(refined|too_fine);
	    else if (cell.flags & too_fine) {
	      if (is_local(cell) && coarsen_cell (point, listc))
		st.nc++;
	      cell.flags &= ~too_fine; // do not coarsen parent
	    }
	  }
	  if (cell.flags & too_fine)
	    cell.flags &= ~too_fine;
	  else if (level > 0 && (aparent(0).flags & too_fine))
	    aparent(0).flags &= ~too_fine;
	  continue;
	}
	else if (is_leaf(cell))
	  continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);
  
  return st;
}
