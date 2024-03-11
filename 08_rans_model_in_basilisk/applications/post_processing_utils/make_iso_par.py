from mpi4py import MPI
from plyfile import PlyData, PlyElement
from ray_triangle import ray_triangle_intersect, Vec3, Ray
import numpy as np
import glob
import csv
from collections import defaultdict
from pathlib import Path
import os
from timeit import default_timer as timer
from scipy.io import savemat


def fill_op(xmem, ymem, dt):
    x = np.frombuffer(xmem)
    y = np.frombuffer(ymem)
    y[~np.isfinite(y)] = x[~np.isfinite(y)]
    y[:] = y


FILL = MPI.Op.Create(fill_op, commute=True)


def readCSV(filename):
    out = defaultdict(list)

    with open(filename, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=' ')
        for line in csv_reader:
            for key, value in line.items():
                out[key].append(value)

    return out


# Utils
def line_distributed_probes(x0, x1, N):
    # xi = [np.linspace(x0[0], x1[0], N), np.linspace(x0[1], x1[1], N), np.linspace(x0[2], x1[2], N)]

    # By default face to - y direction
    dir = Vec3(0, -1, 0)

    rays = [Ray(orig=Vec3(x, y, z), direction=dir) for x, y, z in zip(np.linspace(x0[0], x1[0], N),
                                                                      np.linspace(x0[1], x1[1], N),
                                                                      np.linspace(x0[2], x1[2], N))]
    return rays


def batch_ray_triangle_intersect(rays, faces, vertices, bbox):
    num_rays = len(rays)
    num_faces = faces.shape[0]
    face_index = np.arange(num_faces)

    hitted_face = -1 + np.zeros((num_rays,), dtype=int)
    distance = float("Inf") + np.zeros((num_rays,))
    u = float("nan") + np.zeros((num_rays,))
    v = float("nan") + np.zeros((num_rays,))

    for i in range(num_rays):
        ray = rays[i]

        # early bbox filter (by default we inspect x and z coords)
        face_mark = np.logical_and(
            np.logical_and(bbox['x_min'] <= ray.orig.x, ray.orig.x <= bbox['x_max']),
            np.logical_and(bbox['z_min'] <= ray.orig.z, ray.orig.z <= bbox['z_max']))

        for j in face_index[face_mark]:
            # do the ray casting
            t, ui, vi = ray_triangle_intersect(ray, vertices[faces[j][0]],
                                               vertices[faces[j][1]],
                                               vertices[faces[j][2]])

            if 0 <= t < distance[i]:
                distance[i] = t
                hitted_face[i] = j
                u[i] = ui
                v[i] = vi

    
    return hitted_face, distance, u, v


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("[MPI environment] rank = %d, size = %d.\n" % (rank, size))
infiles = glob.glob("postProcessing/result.*.ply")
optfiles = glob.glob("postProcessing/result.*.csv")
if rank==0:
    print("Found "+ str(len(infiles))+" ply files.")
inputArray = ["u_x", "u_y", "p"]
N = len(infiles)

# Setup rays
y_offset = 0.6
x0 = [3.79, y_offset, 0.0]
x1 = [40.0, y_offset, 0.0]
num_probes = round(40 / 0.01)
rays = line_distributed_probes(x0, x1, num_probes)

# Initial output fields
output_file = "eta_withfields.mat"
corrupted = 0

nbsteps = np.arange(N)
nsub = len(nbsteps[rank::size])
print(nbsteps[rank::size])
fields = {
    "t": np.nan + np.zeros((nsub,)),
    "eta": np.nan + np.zeros((nsub, num_probes)),
}
for arr in inputArray:
    fields[arr] = np.nan + np.zeros((nsub, num_probes))

if rank == 0:
    t_start = timer()

for i, ii in zip(nbsteps[rank::size], range(nsub)):
    # Read ploydata
    plydata = PlyData.read(infiles[i])

    vertices = [Vec3(x, y, z) for x, y, z in
                zip(plydata['vertex']['x'], plydata['vertex']['y'], plydata['vertex']['z'])]
    faces = np.vstack(plydata['face'].data['vertex_index'])
    bbox = {"x_min": np.min(plydata['vertex']['x'][faces], axis=1),
            "y_min": np.min(plydata['vertex']['y'][faces], axis=1),
            "z_min": np.min(plydata['vertex']['z'][faces], axis=1),
            "x_max": np.max(plydata['vertex']['x'][faces], axis=1),
            "y_max": np.max(plydata['vertex']['y'][faces], axis=1),
            "z_max": np.max(plydata['vertex']['z'][faces], axis=1)
            }

    # display some information
    print(str(i) + "/" + str(N))

    fields["t"][ii] = float(infiles[i].rstrip('.ply').split("result.")[-1])

    # Read optional field data
    data = readCSV(infiles[i].rstrip('.ply')+".csv")

    # Perform multi - raycast
    # t_start = timer()
    hitted_faces, distance, w0, w1 = batch_ray_triangle_intersect(rays, faces, vertices, bbox)
    distance -= y_offset
    distance[np.abs(distance) > 1e10] = np.nan
    fields["eta"][ii, :] = distance

    vi0 = faces[hitted_faces][:, 0]
    vi1 = faces[hitted_faces][:, 1]
    vi2 = faces[hitted_faces][:, 2]

    for arr in inputArray:
        fields[arr][ii, :] = w0 * np.array(data[arr], dtype=float)[vi0] + \
                            w1 * np.array(data[arr], dtype=float)[vi1] + \
                            (1.0 - w0 - w1) * np.array(data[arr], dtype=float)[vi2]
        fields[arr][ii, :][np.abs(fields[arr][ii, :]) > 1e10] = np.nan

comm.Barrier()


gafields = comm.gather(fields, root=0)
outfields = {}

comm.Barrier()
if rank==0:
    
    outfields['t'] = np.hstack([gafields[i]['t'] for i in range(len(gafields))])
    outfields['eta'] = np.vstack([gafields[i]['eta'] for i in range(len(gafields))])
    for arr in inputArray:
        outfields[arr] = np.vstack([gafields[i][arr] for i in range(len(gafields))])

    # sort the array
    it = np.argsort(outfields['t'])
    outfields['t'] = outfields['t'][it]
    outfields['eta'] = outfields['eta'][it, :]
    for arr in inputArray:
        outfields[arr] = outfields[arr][it, :]

    # print(outfields['t'].shape,outfields['eta'].shape)
    savemat(output_file, outfields)

    t_end = timer()
    t_total = t_end - t_start

    print('  Total time:\t\t\t  %6.2f seconds' % (t_total))

