void output_vtk2d (scalar * list, int n, FILE * fp, bool linear)
{
  fputs ("# vtk DataFile Version 2.0\n"
	 "Basilisk\n"
	 "ASCII\n"
	 "DATASET STRUCTURED_GRID\n", fp);
  fprintf (fp, "DIMENSIONS %d %d 1\n", n, n);
  fprintf (fp, "POINTS %d double\n", n*n);

  double fn = n;
  double Delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double x = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double y = Delta*j + Y0 + Delta/2.;
      fprintf (fp, "%g %g 0\n", x, y);
    }
  }
  fprintf (fp, "POINT_DATA %d\n", n*n);
  for (scalar s in list) {
    fprintf (fp, "SCALARS %s double\n", s.name);
    fputs ("LOOKUP_TABLE default\n", fp);
    double fn = n;
    double Delta = L0/fn;
    for (int i = 0; i < n; i++) {
      double x = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
	double y = Delta*j + Y0 + Delta/2., v;
	if (linear)
	  v = interpolate (s, x, y);
	else {
	  Point point = locate (x, y);
	  v = point.level >= 0 ? val(s) : nodata;
	}
	fprintf (fp, "%g\n", v);
      }
    }
  }
  fflush (fp);
}

void output_vtk3d (scalar * list, int n, FILE * fp, bool linear)
{
  fputs ("# vtk DataFile Version 2.0\n"
	 "Basilisk\n"
	 "ASCII\n"
	 "DATASET STRUCTURED_GRID\n", fp);
  fprintf (fp, "DIMENSIONS %d %d %d\n", n, n, n);
  fprintf (fp, "POINTS %d double\n", n*n*n);

  double fn = n;
  double Delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double xi = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < n; j++) {
      double yi = Delta*j + Y0 + Delta/2.;
      for (int k = 0; k < n; k++) {
        double zi = Delta*k + Z0 + Delta/2.;
        fprintf (fp, "%g %g %g\n", xi, yi, zi);
      }
    }
  }
  fprintf (fp, "POINT_DATA %d\n", n*n*n);
  for (scalar s in list) {
    fprintf (fp, "SCALARS %s double\n", s.name);
    fputs ("LOOKUP_TABLE default\n", fp);
    double fn = n;
    double Delta = L0/fn;
    for (int i = 0; i < n; i++) {
      double xi = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < n; j++) {
	double yi = Delta*j + Y0 + Delta/2.;
        for (int k = 0; k < n; k++) {
          double zi = Delta*k + Z0 + Delta/2., v;
	  Point point = locate (xi, yi, zi);
	  v = point.level >= 0 ? val(s) : nodata;
	  fprintf (fp, "%g\n", v);
        }
      }
    }
  }
  fflush (fp);
}
