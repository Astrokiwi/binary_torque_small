from gtools import gizmo_tools
import numpy as np
import tempfile
import os
import h5py
import pynbody.plot.sph as sph
import pynbody.plot.generic as generic
from pynbody import filt

# Parameters for dumps - times (t_orbit)
runs = [["norad_prod", ["small_ecc","small_circ"]],
        ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"]],
        ["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"]],
        ["norad_unary", ["close"]],
        ["rad_unary", ["close"]]
        ]
t_orbit = 13.7e-3 # Myr
# t_images = [[0,2,5,10,50,100],
#             [0,2,5,10],
#             [0,2,5,10]]

t_images = [[2,20,100],
            [1,2,5],
            [1,2,5],
            [2, 20, 100],
            [1,2,5],
            ]

widths = [0.2,4.]

res = 256

if __name__ == '__main__':
    # Get snapshots for times
    isnap_images = \
            [
                [np.searchsorted(
                    np.array(gizmo_tools.snapshot_times(run_dir,run_name))/t_orbit,
                    t_image)
                for run_name in run_names]
            for t_image,(run_dir, run_names) in zip(t_images,runs)]

    # Render & save plots
    f_out = h5py.File("../analysis_out/image_maps_evolution.hdf5", "w")

    for isnap_runs,(run_dir, run_names) in zip(isnap_images,runs):
        gizmoDir = gizmo_tools.getGizmoDir(run_dir)
        with tempfile.NamedTemporaryFile() as tf :
            for run_name,snap_list in zip(run_names,isnap_runs):
                fullDir = os.path.join(gizmoDir, run_dir, run_name)
                for ii,isnap in enumerate(snap_list):
                    print(f"Rendering {run_name} {isnap}")
                    infile = f"{fullDir}/snapshot_{isnap:03d}.hdf5"
                    header,snap = gizmo_tools.load_gizmo_nbody(run_dir,run_name,f"{isnap:03d}")
                    for iwidth,width in enumerate(widths):
                        m = sph.image(snap.gas, qty="rho", units="g cm^-3", width=f"{width} pc", resolution=res, noplot=True, approximate_fast=True)
                        key = f"{run_dir}_{run_name}_{ii}_{iwidth}_face"
                        f_out[key] = m # g/cm^2
                        # centre = snap[filt.disc(f"{width*2} pc",f"{width*2} pc")]
                        snap.rotate_x(90)
                        m = sph.image(snap.gas, qty="rho", units="g cm^-3", width=f"{width} pc", resolution=res, noplot=True, approximate_fast=True)
                        # m = generic.hist2d(centre["rxy"],centre["z"],gridsize=(res,res),weights=centre["mass"],make_plot=False)
                        # cell size is
                        snap.rotate_x(270)
                        key = f"{run_dir}_{run_name}_{ii}_{iwidth}_edge"
                        f_out[key] = m
    f_out.close()