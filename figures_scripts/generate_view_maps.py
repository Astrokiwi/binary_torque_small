from gtools import gizmo_tools
import numpy as np
from multiprocessing import Pool
from gtools import frame
import os
import tqdm
import h5py

# Parameters for dumps - times (t_orbit)
# runs = [["norad_prod", ["small_ecc","small_circ"], ["norad_ecc","norad_circ"]],
#         ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"], ["rad_ecc_latecont","rad_circ_latecont"]],
#         ["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"], ["rad_ecc_earlycont","rad_circ_earlycont"]],
#         ["norad_unary", ["close"],["norad_solo"]],
#         ["rad_unary", ["close"],["rad_solo"]]
#         ]
runs = [
        ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"], ["rad_ecc_latecont","rad_circ_latecont"]],
        ["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"], ["rad_ecc_earlycont","rad_circ_earlycont"]],
        ["rad_unary", ["close"],["rad_solo"]]
        ]
t_orbit = 13.7e-3 # Myr
# t_images = [[0,2,5,10,50,100],
#             [0,2,5,10],
#             [0,2,5,10]]

# t_images = [[2,20,100],
#             [1,2,5],
#             [1,2,5],
#             [2, 10, 18],
#             [1,3]
#             ]
t_images = [[2],
            [2],
            [2]
            ]

widths = [0.2,1.]

res = 128
# res = 64

major_axis = 0.035

rots = np.linspace(0.,np.pi/2.,12) # representative

nprocs = 32

def gen_image(image_config):
    run_code, run_dir, run_name, isnap, width, rot = image_config
    # print(run_code)
    gizmoDir = gizmo_tools.getGizmoDir(run_dir)
    fullDir = os.path.join(gizmoDir, run_dir, run_name)
    infile = f"{fullDir}/snapshot_{isnap:03d}.hdf5"
    outfile = f"../quick_out/rot_plots/{run_code}.png"


    retmap = frame.makesph_trhoz_frame(infile, outfile, cmap='plasma', flat=True, ring=False,
                          plot=["view"], L=res, pixsize=4, views=['face'],gaussian=0.16,
                          rot=[0., rot], scale=width, visibleAxes=True,
                          return_maps=True)
    # print("shape=",retmap[0].shape)
    return run_code,retmap

if __name__ == '__main__':
    print("Calculating snapshot times")
    # Get snapshots for times
    isnap_images = \
            [
                [np.searchsorted(
                    np.array(gizmo_tools.snapshot_times(run_dir,run_name))/t_orbit,
                    t_image)
                for run_name in run_names]
            for t_image,(run_dir, run_names, run_realnames) in zip(t_images,runs)]

    # Render & save plots

    # build list of images
    print("Building list")
    image_list = []
    for isnap_runs,(run_dir, run_names,run_realnames) in zip(isnap_images,runs):
        for run_name,snap_list in zip(run_names,isnap_runs):
            for ii,isnap in enumerate(snap_list):
                for iwidth,width in enumerate(widths):
                    for irot,rot in enumerate(rots):
                        image_list.append([f"{run_dir}_{run_name}_{ii}_{iwidth}_{irot}",run_dir,run_name,isnap,width,rot])

    with Pool(processes=nprocs) as pool:
        retmaps=[]
        for _ in tqdm.tqdm(pool.imap(gen_image,image_list),total=len(image_list)):
            retmaps.append(_)

    print(len(retmaps),"saving")
    f = h5py.File("../analysis_out/image_maps_revolution.hdf5", "w")
    for map_code,image_map in retmaps:
        f[map_code]=image_map
    f.close()
