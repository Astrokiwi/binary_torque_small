from gtools import gizmo_tools
import h5py

import numpy as np
import pynbody

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

print("Loading")
header,snap = gizmo_tools.load_gizmo_nbody("rad_prod","rad_small_circ_earlier","050",load_binary_headers=True,load_vals=["bh_pos"])
cent = snap[pynbody.filt.Sphere("0.2 pc")]

fin = h5py.File("../analysis_out/tau_test.hdf5",'r')
d = np.array(fin["tau_sim"])
fin.close()

fig,sp = plt.subplots(2,1,figsize=(6,12))
h,xedges,yedges,im = sp[0].hist2d(d,cent["AGNDepth"],bins=(200,200))
sp[0].plot(xedges,xedges)

disc_slice = np.abs(cent["z"].in_units("pc"))<1.e-3
disc = cent[disc_slice]
d_disc = d[disc_slice]
img = sp[1].scatter(disc["x"],disc["y"],c=d_disc,s=0.1,marker='.')
fig.colorbar(img,ax=sp[1])

fig.savefig("../figures_out/resample_ray_testviz.png",dpi=100)