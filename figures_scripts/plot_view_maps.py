import numpy as np
import h5py
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from gtools import gizmo_tools

# Parameters for dumps - times (t_orbit)

runs = [
        ["rad_prod","rad_small_ecc_full","rad_ecc_latecont"],
        ["rad_prod","rad_small_circ_full","rad_circ_latecont"],
        ["rad_prod","rad_small_ecc_earlier","rad_ecc_earlycont"],
        ["rad_prod","rad_small_circ_earlier","rad_circ_earlycont"],
        ["rad_unary","close","rad_solo"]
        ]

t_orbit = 13.7e-3 # Myr

widths = [0.2,1.]

res = 128

major_axis = 0.035

widths = np.array(widths)/major_axis

rots = np.degrees(np.linspace(0.,np.pi/2.,12))

# rot_selected = list(range(12))
rot_selected = [0,5,8]

nruns = len(runs)
nwidths = len(widths)
nrots = len(rot_selected)

if __name__ == '__main__':
    code_base = "{run_dir}_{run_name}_{ii}_{iwidth}_{irot}"

    f = h5py.File("../analysis_out/image_maps_revolution.hdf5", "r")

    nx = nruns
    ny = nrots*2
    # figs = [None,None]
    # sps = [None,None]

    scale = 2.5
    cmap = 'inferno'
    box_corner = [-2.2,-2.2]
    box_dims = [4.4,4.4]
    box_color='green'

    # figs[0],sps[0] = plt.subplots(ny,nx,figsize=(scale*nx,scale*ny),sharex=True,sharey=True,constrained_layout=True)
    # figs[1],sps[1] = plt.subplots(ny,nx,figsize=(scale*nx,scale*ny),sharex=True,sharey=True,constrained_layout=True)

    fig,sp = plt.subplots(ny,nx,figsize=(scale*nx,scale*ny),sharex='row',sharey='row',constrained_layout=True)

    for irun,(run_dir,run_name,run_title) in enumerate(runs):
        for iwidth,width in enumerate(widths):
            for irot,rot_index in enumerate(rot_selected):
                ix = irun
                # fig = figs[iwidth]
                # sp = sps[iwidth]
                # iy = irot+nrots*iwidth
                iy = irot*2+iwidth

                ax = sp[iy,ix]
                if iy==0:
                    ax.set_title(run_title)
                if ix==nx-1:
                    ax.yaxis.set_label_position("right")
                    ax.set_ylabel(fr"${int(rots[rot_index]):2d}^\circ$")

                im_map = f[code_base.format(run_dir=run_dir,run_name=run_name,ii=0,iwidth=iwidth,irot=rot_index)][0,:,:]
                im_map = np.flipud(im_map)
                # print(np.min(im_map),np.max(im_map))
                ax.imshow(im_map,
                    extent = [-width / 2, width / 2, -width / 2, width / 2],
                                 cmap=cmap,
                                 vmin=0.,vmax=6.)

    for irun in range(nruns) :
        for irot in range(nrots):
            ix = irun
            iy = irot * 2
            gizmo_tools.box_connected_two_axes(sp[iy, ix], sp[iy+1, ix], box_corner, box_dims,color=box_color,corners=[[0,0],[1,0]])

    f.close()
    fig.savefig("../figures_out/rot_maps.png")
    # figs[0].savefig("../figures_out/rot_maps_close.png")
    # figs[1].savefig("../figures_out/rot_maps_far.png")