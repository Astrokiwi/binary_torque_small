
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0 :
    print("Importing (MPI already imported)")

from gtools import gizmo_tools
import ray_py
import h5py

import numpy as np
import pynbody



if rank == 0 :
    print("Loading")
    header, snap = gizmo_tools.load_gizmo_nbody("rad_prod", "rad_small_circ_earlier", "050", load_binary_headers=True,
                                                load_vals=["bh_pos"])
    cent = snap#[pynbody.filt.Sphere("0.1 pc")]
    print(len(cent))
    # plane = cent[pynbody.filt.Disc("0.1 pc","0.001 pc")]
    r = np.array(cent["pos"])
    smooth = np.array(cent["smooth"])
    opac = np.array(cent["AGNOpacity"])
    r_agn = 1.e-3*np.array(header["BH_pos_1"])
else:
    r_agn = None
    r=None
    smooth=None
    opac=None

r_agn = comm.bcast(r_agn, root=0)

if rank == 0 :
    print("Raytracing")

d = ray_py.calc_depths(r,smooth,opac,r_agn=r_agn,mass=1.e-14)

if rank == 0 :
    print("Dumping",len(d))
    fout = h5py.File("../analysis_out/tau_test.hdf5",'w')
    fout["tau_sim"] = d
    fout.close()
