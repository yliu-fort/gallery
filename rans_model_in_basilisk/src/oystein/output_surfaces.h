/**
Various output functions for extracting a surface from field data using the build in VOF PLIC surface or the isosurface reconstruction in bview
Three different formats are supported:
- .ply
- .vtu
- .vtp
all of which can be loaded into Paraview. for the two later alternatives, field data may be added by using the functions with "_w_fielddata".
This is an updated version of the previous fractions_output.h

Author: Oystein Lande
Date: July 2019

Modified by: Yuxuan Liu
Date: Jan 2022
*/

#include "geometry.h"
#include "fractions.h"

#if dimension == 1
coord mycs (Point point, scalar c) {
    coord n = {1.};
    return n;
}
#elif dimension == 2
# include "myc2d.h"
#else // dimension == 3
# include "myc.h"
#endif


/**
 This function outputs the VOF surface to a .ply file
 */
struct OutputFacets_scalar {
    scalar c;
    FILE * fp;     // optional: default is stdout
    FILE * fp2;     // optional: default is stderr
    scalar * list;  // List of scalar fields to include when writing vtu surface to file
    vector * vlist; // List of vector fields to include.
    face vector s; // optional: default is none
};

int replacechar(char *str, char orig, char rep) {
    char *ix = str;
    int n = 0;
    while((ix = strchr(ix, orig)) != NULL) {
        *ix++ = rep;
        n++;
    }
    return n;
}

void output_ply_mpi_compatible (struct OutputFacets_scalar p)
{

    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;
    if (!p.fp) p.fp = stdout;
    if (!p.fp2) p.fp2 = stderr;

    // print header text
    fputs ("ply\n", p.fp);
    fputs ("format ascii 1.0\n", p.fp);

	int total_nverts = 0;
	int total_nfacets = 0;
	long offset = 0;
	long offset2 = 0;

    int world_size = 1;
#if _MPI
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Barrier(MPI_COMM_WORLD);
#endif

// start from processor 0, write the field by sequence
// once completed, broadcast the total verts and facets to other processors

	int* nverts = calloc(world_size, sizeof(int));
	int* nfacets = calloc(world_size, sizeof(int));

    foreach(reduction(+:nverts[:world_size]) reduction(+:nfacets[:world_size]))
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs (point, c);
            else {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension() {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert (nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha (c[], n);

            coord v[12];
            int ntris = 0;
#if dimension == 3
            int m = facets (n, alpha, v, 1.);
            ntris = m - 2;
            m = ntris * 3;
#elif dimension == 2
            // in 2d, two intersection points are defined which is v[0](v0.x,v0.y,0)
            // and v[1] (v1.x,v1.y,0)
            // first project to 4 vertices to 3d space:
            // v[0] = (v0.x,v0.y,-z0) v[3] = (v0.x,v0.y,+z0)
            // v[1] = (v1.x,v1.y,-z0) v[2] = (v1.x,v1.y,+z0)
            // connect to 2 triangles in z-direction defined by
            // (v0, v1,v2), (v0, v2, v3)
            int m = facets (n, alpha, v);
            if (m == 2) {
            	double v0x = v[0].x;
            	double v0y = v[0].y;
            	double v1x = v[1].x;
            	double v1y = v[1].y;

            	v[0].x = v0x;
            	v[0].y = v0y;
            	v[0].z =-0.5;

            	v[1].x = v1x;
            	v[1].y = v1y;
            	v[1].z =-0.5;

            	v[2].x = v1x;
            	v[2].y = v1y;
            	v[2].z = 0.5;

            	v[3].x = v0x;
            	v[3].y = v0y;
            	v[3].z = 0.5;

            	m = 4;
            	ntris = m - 2;
            	m = ntris * 3;
            }
            else{
            	m = 0;
            	ntris = 0;
            }

#endif
            for (int i = 0; i < m; i++) {
                nverts[pid()] ++;
            }
            if (ntris > 0) {
                nfacets[pid()] += ntris;
            }
        }

#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&(nverts[pid()]), &total_nverts, 1, MPI_INT, MPI_SUM, 0,
	           MPI_COMM_WORLD);
	MPI_Reduce(&(nfacets[pid()]), &total_nfacets, 1, MPI_INT, MPI_SUM, 0,
	           MPI_COMM_WORLD);
#else
	total_nverts = *nverts;
	total_nfacets = *nfacets;
#endif

	// pid 0 writes the headers
    if (pid() == 0){
        fprintf (p.fp, "element vertex %i\n", total_nverts);
        fputs ("property float x\n", p.fp);
        fputs ("property float y\n", p.fp);
        fputs ("property float z\n", p.fp);

        // Write headers of interpolated fields to fp2
        for (scalar s in p.list) {
            fprintf (p.fp2,"%s ", s.name);
        }
        for (vector v in p.vlist) {
        	char vname[80];
        	strcpy(vname, v.x.name);
        	replacechar(vname, '.','_');
            fprintf (p.fp2,"%s ", vname);
            strcpy(vname, v.y.name);
			replacechar(vname, '.','_');
            fprintf (p.fp2,"%s ", vname);
#if dimension == 3
            strcpy(vname, v.z.name);
			replacechar(vname, '.','_');
            fprintf (p.fp2,"%s ", vname);
#endif
        }
        fprintf (p.fp2,"\n");

        fprintf (p.fp, "element face %i\n", total_nfacets);
        fputs ("property list uchar int vertex_index\n", p.fp);
        fputs ("end_header\n", p.fp);

        fflush (p.fp);
        fflush (p.fp2);
        offset = ftell(p.fp);
        offset2 = ftell(p.fp2);
    }

#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	//MPI_Bcast(&offset, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	if (pid() != 0) {
	    MPI_Recv(&offset, 1, MPI_LONG, pid() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&offset2, 1, MPI_LONG, pid() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	fseek (p.fp, offset, SEEK_SET);
	fseek (p.fp2, offset2, SEEK_SET);
#endif
        //int facet_num[nfacets]; // Since we triangulized the facets
        //int* facet_num = malloc(nfacets*sizeof(*facet_num));

        //int ifacet = 0;
        int ivert = 0;
	    //int* ifacet = calloc(world_size, sizeof(int));
	    //int* ivert = calloc(world_size, sizeof(int));

        foreach()
            if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                coord n;
                if (!s.x.i)
                    // compute normal from volume fraction
                    n = mycs (point, c);
                else {
                    // compute normal from face fractions
                    double nn = 0.;
                    foreach_dimension() {
                        n.x = s.x[] - s.x[1];
                        nn += fabs(n.x);
                    }
                    assert (nn > 0.);
                    foreach_dimension()
                        n.x /= nn;
                }
                double alpha = plane_alpha (c[], n);

                coord v[12];
                int ntris = 0;
#if dimension == 3
				int m = facets (n, alpha, v, 1.);
				ntris = m - 2;
				m = ntris * 3;
#elif dimension == 2
				// in 2d, two intersection points are defined which is v[0](v0.x,v0.y,0)
				// and v[1] (v1.x,v1.y,0)
				// first project to 4 vertices to 3d space:
				// v[0] = (v0.x,v0.y,-z0) v[3] = (v0.x,v0.y,+z0)
				// v[1] = (v1.x,v1.y,-z0) v[2] = (v1.x,v1.y,+z0)
				// connect to 2 triangles in z-direction defined by
				// (v0, v1,v2), (v0, v2, v3)
                int vz[4];
				int m = facets (n, alpha, v);
				if (m == 2) {
					double v0x = v[0].x;
					double v0y = v[0].y;
					double v1x = v[1].x;
					double v1y = v[1].y;

					v[0].x = v0x;
					v[0].y = v0y;
					v[0].z =-0.5;

					v[1].x = v1x;
					v[1].y = v1y;
					v[1].z =-0.5;

					v[2].x = v1x;
					v[2].y = v1y;
					v[2].z = 0.5;

					v[3].x = v0x;
					v[3].y = v0y;
					v[3].z = 0.5;

					m = 4;
					ntris = m - 2;
					m = ntris * 3;
				}
				else{
	            	m = 0;
	            	ntris = 0;
	            }

#endif

            	for (int j = 0; j < ntris; j++) {

					double v0, v1, v2;

					for (int k = 0; k < 3; k++) {
						int ind = k == 0?0:j+k;
						v0 = x + v[ind].x*Delta;
						v1 = y + v[ind].y*Delta;
						v2 = z + v[ind].z*Delta;
						fprintf (p.fp, "%g %g %g\n", v0, v1, v2);

						// Write interpolated fields to fp2
                        
                        Point point = locate (v0, v1, v2);
                        
                        if (point.level >= 0) {
                            for (scalar s in p.list)
                                fprintf (p.fp2,"%g ", interpolate_linear (point, (struct _interpolate){s, v0, v1, v2}));
                            for (vector v in p.vlist)
                                foreach_dimension(){
                                    fprintf (p.fp2,"%g ", interpolate_linear (point, (struct _interpolate){v.x, v0, v1, v2}));}
                        }
                        else {
                            for (scalar s in p.list)
                                fprintf (p.fp2,"%g ", nodata);
                            for (vector v in p.vlist)
                                foreach_dimension(){
                                    fprintf (p.fp2,"%g ", nodata);}
                        }
				        fprintf (p.fp2, "\n");
					}

					//if (m > 0) {
							//facet_num[ifacet] = 3;
					//		ifacet[pid()] ++;
					//}
                }
            }
        fflush (p.fp);
        fflush (p.fp2);
        offset = ftell(p.fp);
        offset2 = ftell(p.fp2);
#if _MPI
	MPI_Send(&offset, 1, MPI_LONG, (pid() + 1) % world_size, 0, MPI_COMM_WORLD);
	MPI_Send(&offset2, 1, MPI_LONG, (pid() + 1) % world_size, 0, MPI_COMM_WORLD);
	if (pid() == 0) {
	    MPI_Recv(&offset, 1, MPI_LONG, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&offset2, 1, MPI_LONG, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	int vert_offset = 0;

#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if (pid() != 0) {
	    MPI_Recv(&vert_offset, 1, MPI_INT, pid() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&offset, 1, MPI_LONG, pid() - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	fseek (p.fp, offset, SEEK_SET);
#endif
            // print face list
        	// num_vertex v0 v1 v2 ...
            for (int ifacet = 0; ifacet < nfacets[pid()]; ifacet++) {
                fprintf (p.fp, "%i ", 3);
                for (int iv = 0; iv < 3; iv ++) {
                    fprintf (p.fp, "%i ", vert_offset+ivert);
                    ivert ++;
                }
                fputc ('\n', p.fp);
            }
            fflush (p.fp);
            offset = ftell(p.fp);
#if _MPI
	vert_offset += ivert;
	MPI_Send(&vert_offset, 1, MPI_INT, (pid() + 1) % world_size, 0, MPI_COMM_WORLD);
	MPI_Send(&offset, 1, MPI_LONG, (pid() + 1) % world_size, 0, MPI_COMM_WORLD);
	if (pid() == 0) {
	    MPI_Recv(&vert_offset, 1, MPI_INT, world_size - 1, 0,
	             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(&offset, 1, MPI_LONG, world_size - 1, 0,
	             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif

            fflush (p.fp);
            fflush (p.fp2);
            #if defined(_OPENMP)
            omp_set_num_threads(num_omp);
            #endif

            free(nverts);
            free(nfacets);
}

#include "ray_math.h"
struct OutputFacets_raycasting {
    scalar c;
    Ray * rays; // rays to cast
    size_t num_rays; // number of rays

    FILE * fp;     // optional: default is stdout
    scalar * list;  // List of scalar fields to include when writing vtu surface to file
    vector * vlist; // List of vector fields to include.
    face vector s; // optional: default is none
};

void reconstruct_interface (struct OutputFacets_raycasting p)
{

    #if defined(_OPENMP)
    int num_omp = omp_get_max_threads();
    omp_set_num_threads(1);
    #endif
    scalar c = p.c;
    face vector s = p.s;

#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

// start from processor 0, write the field by sequence
// once completed, broadcast the total verts and facets to other processors
    int num_fields = 1; // f
    for (scalar s in p.list)
        num_fields++;
#if dimension == 2
    for (vector v in p.vlist){num_fields+=2;}
#else
    for (vector v in p.vlist){num_fields+=3;}
#endif

    double **result = calloc(p.num_rays, sizeof(double*));
    for(int i=0; i < p.num_rays; i++)
    {
        result[i] = calloc(num_fields, sizeof(double));
        for(int j=0; j < num_fields; j++)
            result[i][j] = nodata;
    }

	int nverts = 0;
	int nfacets = 0;

    foreach(noauto)
        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
            coord n;
            if (!s.x.i)
                // compute normal from volume fraction
                n = mycs (point, c);
            else {
                // compute normal from face fractions
                double nn = 0.;
                foreach_dimension() {
                    n.x = s.x[] - s.x[1];
                    nn += fabs(n.x);
                }
                assert (nn > 0.);
                foreach_dimension()
                    n.x /= nn;
            }
            double alpha = plane_alpha (c[], n);

            coord v[12];
            int ntris = 0;
#if dimension == 3
            int m = facets (n, alpha, v, 1.);
            ntris = m - 2;
            m = ntris * 3;
#elif dimension == 2
            // in 2d, two intersection points are defined which is v[0](v0.x,v0.y,0)
            // and v[1] (v1.x,v1.y,0)
            // first project to 4 vertices to 3d space:
            // v[0] = (v0.x,v0.y,-z0) v[3] = (v0.x,v0.y,+z0)
            // v[1] = (v1.x,v1.y,-z0) v[2] = (v1.x,v1.y,+z0)
            // connect to 2 triangles in z-direction defined by
            // (v0, v1,v2), (v0, v2, v3)
            int m = facets (n, alpha, v);
            if (m == 2) {
            	double v0x = v[0].x;
            	double v0y = v[0].y;
            	double v1x = v[1].x;
            	double v1y = v[1].y;

            	v[0].x = v0x;
            	v[0].y = v0y;
            	v[0].z =-0.5;

            	v[1].x = v1x;
            	v[1].y = v1y;
            	v[1].z =-0.5;

            	v[2].x = v1x;
            	v[2].y = v1y;
            	v[2].z = 0.5;

            	v[3].x = v0x;
            	v[3].y = v0y;
            	v[3].z = 0.5;

            	m = 4;
            	ntris = m - 2;
            	m = ntris * 3;
            }
            else{
            	m = 0;
            	ntris = 0;
            }

#endif
            for (int i = 0; i < m; i++) {
                nverts ++;
            }
            if (ntris > 0) {
                nfacets += ntris;
            }
        }

#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

        foreach(noauto)
            if (c[] > 1e-6 && c[] < 1. - 1e-6) {
                coord n;
                if (!s.x.i)
                    // compute normal from volume fraction
                    n = mycs (point, c);
                else {
                    // compute normal from face fractions
                    double nn = 0.;
                    foreach_dimension() {
                        n.x = s.x[] - s.x[1];
                        nn += fabs(n.x);
                    }
                    assert (nn > 0.);
                    foreach_dimension()
                        n.x /= nn;
                }
                double alpha = plane_alpha (c[], n);

                coord v[12];
                int ntris = 0;
#if dimension == 3
				int m = facets (n, alpha, v, 1.);
				ntris = m - 2;
				m = ntris * 3;
#elif dimension == 2
				// in 2d, two intersection points are defined which is v[0](v0.x,v0.y,0)
				// and v[1] (v1.x,v1.y,0)
				// first project to 4 vertices to 3d space:
				// v[0] = (v0.x,v0.y,-z0) v[3] = (v0.x,v0.y,+z0)
				// v[1] = (v1.x,v1.y,-z0) v[2] = (v1.x,v1.y,+z0)
				// connect to 2 triangles in z-direction defined by
				// (v0, v1,v2), (v0, v2, v3)
                int vz[4];
				int m = facets (n, alpha, v);
				if (m == 2) {
					double v0x = v[0].x;
					double v0y = v[0].y;
					double v1x = v[1].x;
					double v1y = v[1].y;

					v[0].x = v0x;
					v[0].y = v0y;
					v[0].z =-0.5;

					v[1].x = v1x;
					v[1].y = v1y;
					v[1].z =-0.5;

					v[2].x = v1x;
					v[2].y = v1y;
					v[2].z = 0.5;

					v[3].x = v0x;
					v[3].y = v0y;
					v[3].z = 0.5;

					m = 4;
					ntris = m - 2;
					m = ntris * 3;
				}
				else{
	            	m = 0;
	            	ntris = 0;
	            }

#endif

            	for (int j = 0; j < ntris; j++) {

					Vec3 _v[3];
                    float x_min= 1e30;
                    float x_max=-1e30;
                    float z_min= 1e30;
                    float z_max=-1e30;
                    // print vertices of the triangle
					for (int k = 0; k < 3; k++) {
						int ind = k == 0?0:j+k;
						_v[k].x = x + v[ind].x*Delta;
						_v[k].y = y + v[ind].y*Delta;
						_v[k].z = z + v[ind].z*Delta;
                        x_min = x_min < _v[k].x ? x_min : _v[k].x;
                        x_max = x_max > _v[k].x ? x_max : _v[k].x;
                        z_min = z_min < _v[k].z ? z_min : _v[k].z;
                        z_max = z_max > _v[k].z ? z_max : _v[k].z;
					}

                    // loop all the rays to perform ray casting
                    for(int ir=0; ir < p.num_rays; ir++)
                    {
                        // Pre-filter for rays fall outside of bbox
                        if(p.rays[ir].orig.x < x_min || p.rays[ir].orig.x > x_max ||
                           p.rays[ir].orig.z < z_min || p.rays[ir].orig.z > z_max)
                            continue;

                        Hit hit = ray_triangle_intersect(p.rays[ir], _v[0], _v[1], _v[2]);
                        if(hit.hit && hit.distance >=0)
                        {
                            // Write the surface elevation
                            double elevation_from_zero_level = p.rays[ir].orig.y - hit.distance;

                            // if the surface elevation is written then return.
                            if (result[ir][0] != nodata && result[ir][0] > elevation_from_zero_level)
                                continue;

                            result[ir][0] = elevation_from_zero_level;

                            // interpolate fields
                            double w0 = hit.u;
                            double w1 = hit.v;
                            double hit_0 = w0 * _v[0].x + w1 * _v[1].x + (1.0 - w0 - w1) * _v[2].x;
                            double hit_1 = w0 * _v[0].y + w1 * _v[1].y + (1.0 - w0 - w1) * _v[2].y;
                            double hit_2 = w0 * _v[0].z + w1 * _v[1].z + (1.0 - w0 - w1) * _v[2].z;
                            Point point = locate (hit_0, hit_1, hit_2);
                            int iff = 1;
                            if (point.level >= 0) {
                                for (scalar s in p.list)
                                    result[ir][iff++] = interpolate_linear (point, (struct _interpolate){s, hit_0, hit_1, hit_2});
                                for (vector v in p.vlist)
                                    foreach_dimension(){
                                        result[ir][iff++] = interpolate_linear (point, (struct _interpolate){v.x, hit_0, hit_1, hit_2});
                                    }
                            }
                        }
                    }
                }
            }

#if _MPI
	MPI_Barrier(MPI_COMM_WORLD);
    // Merge results
    for(int ir=0; ir < p.num_rays; ir++){
        if(pid()==0){
            MPI_Reduce(MPI_IN_PLACE, result[ir], num_fields, MPI_DOUBLE, MPI_MIN, 0,
                        MPI_COMM_WORLD);
        }else{
            MPI_Reduce(result[ir], result[ir], num_fields, MPI_DOUBLE, MPI_MIN, 0,
                        MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if(pid() == 0){
        // Write the header
        fprintf (p.fp,"x ");
#if dimension == 3
        fprintf(p.fp,"z ");
#endif
        fprintf (p.fp,"eta ");
        for (scalar s in p.list) {
            fprintf (p.fp,"%s ", s.name);
        }
        for (vector v in p.vlist) {
        	char vname[80];
        	strcpy(vname, v.x.name);
        	replacechar(vname, '.','_');
            fprintf (p.fp,"%s ", vname);
            strcpy(vname, v.y.name);
			replacechar(vname, '.','_');
            fprintf (p.fp,"%s ", vname);
#if dimension == 3
            strcpy(vname, v.z.name);
			replacechar(vname, '.','_');
            fprintf (p.fp,"%s ", vname);
#endif
        }
        fprintf (p.fp,"\n");

    // Write the results
        for(int i=0; i < p.num_rays; i++)
        {
            fprintf (p.fp, "%g ", p.rays[i].orig.x);
#if dimension == 3
            fprintf (p.fp, "%g ", p.rays[i].orig.z);
#endif
            for(int j = 0;j < num_fields; j++)
                fprintf (p.fp, "%g ", result[i][j] == nodata ? NAN : result[i][j]);
            fprintf (p.fp, "\n");
        }
}
for(int i=0; i < p.num_rays; i++)
    free(result[i]);
free(result);

#if defined(_OPENMP)
    omp_set_num_threads(num_omp);
#endif
}