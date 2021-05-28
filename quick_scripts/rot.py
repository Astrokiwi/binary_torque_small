print("Importing")

# Import libraries to do our magic
import numpy as np
import sys
import os
from multiprocessing import Pool

from gtools import frame

from gtools import gizmo_tools
import tqdm

print("Running")

# run_id = sys.argv[1]
# output_dir = sys.argv[2]
#
# snapx = int(sys.argv[3])


nprocs = 128




rots = np.linspace(0.,np.pi/2.,24) # representative
# rots = np.linspace(0.,np.pi*2.,360.) # smooth spin
# rots = np.linspace(0.,np.pi,180.) # smooth spin, no-loop
# rots = np.linspace(0.,np.pi,64) # smooth-ish, no-loop
# rots = np.linspace(0.,np.pi*2.,2) # test

titlesuffixes = [r" $\phi=%2d^\circ$"%(theta*180./np.pi) for theta in rots]

pic_dir = gizmo_tools.getPicDir()
movieDir = gizmo_tools.getMovieDir()


for run_id,output_dir,snapx in [
        ["rad_unary","close",20],
        ["rad_prod", "rad_small_circ_earlier", 20],
        ["rad_prod", "rad_small_circ", 20],
        ["rad_prod", "rad_small_ecc_earlier", 20],
        ["rad_prod", "rad_small_ecc", 20]
    ]:

    gizmoDir = gizmo_tools.getGizmoDir(run_id)
    fullDir = gizmoDir+"/"+run_id+"/"+output_dir
    outfiles = [pic_dir+"/sphrotplot"+run_id+output_dir+"%03d_%03d.png"%(snapx,irot) for irot in range(len(rots))]
    os.system("rm "+pic_dir+"/sphrotplot"+run_id+output_dir+"%03d"%snapx+"_???.png")
    infile = fullDir+"/snapshot_%03d.hdf5" % snapx


    def frame_i(irot):
        return frame.makesph_trhoz_frame(infile, outfiles[irot], cmap='plasma', flat=True, ring=False,
                                     plot=["view"], L=200, pixsize=4, gaussian=0.16, views=['face'],
                                     rot=[0., rots[irot]], scale=.5, visibleAxes=True)

    with Pool(processes=nprocs) as pool:
        maps=[]
        for _ in tqdm.tqdm(pool.imap_unordered(frame_i,range(len(rots))),total=len(rots)):
            maps.append(_)

    for result in maps:
        if isinstance(result, frame.ExceptionWrapper):
            result.re_raise()

    cmd = "ffmpeg -y -r 24 -i "+pic_dir+"/sphrotplot"+run_id+output_dir+"%03d_"%snapx+"%03d.png -c:v mpeg4 -q:v 1 "+movieDir+"/rotateview_"+run_id+"_"+output_dir+"_%03d"%snapx+".mp4"

    os.system(cmd)
