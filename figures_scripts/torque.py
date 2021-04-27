import h5py
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File("../analysis_out/tidy_torque.hdf5",'r')
fig,sp = plt.subplots(2,2,figsize=(9,6),sharex=True,sharey=False,constrained_layout=True)

plot_prams = {
    "acc": {
        "title" : r"$\Delta L_{acc}/L_0$"
        },
    "grav" : {
        "title" : r"$\Delta L_{grav}/L_0$"
        }

    }

for isim in range(2):
    for bh in ["1", "2"] :
        for irun in range(3) :
            time = np.array(f[f"time_{irun}_{isim}"])
            for iy,(plot_key,ax_pram) in enumerate(plot_prams.items()):
                ax = sp[iy,isim]
                h5_key = f"BH_{plot_key}_J_{bh}_{irun}_{isim}"
                y = -np.array(f[h5_key]) # minus sign so that angular momentum of binary is positive
                ax.plot(time,y)

for isp,ax_pram in enumerate(plot_prams.values()):
    sp[isp,0].set_ylabel(ax_pram["title"])

for ix in range(2):
    sp[-1,ix].set_xlabel('t (orbits)')

sp[0,0].set_title("ecc")
sp[0,1].set_title("circ")

fig.savefig("../figures_out/torque.pdf")
plt.close(fig)
f.close()