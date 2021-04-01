import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from astropy.io import fits
import numpy as np
import os

simulations = ["restest"]
# simulations = ["binary","unary"]
# imagers = ["K","L","12"]
# imagers = ["K","bands"]
imagers = ["bands"]
# angles = ["face","nearface","edge","nearedge"]
# angles = [f"{inc:02d}" for inc in np.linspace(0,90,11).astype(np.int)]
angles = ["00","45","90"]

nx = len(imagers)*len(simulations)
ny = len(angles)

with open("../local_config/skirt_data_dir") as f:
    skirt_data_dir = f.readline()

# file_base = "{skirt_data_dir}/{sim}_lowres_{imgr}_{angle}_total.fits"
file_base = "{skirt_data_dir}/test_data/big_{sim}_{imgr}_{angle}_total.fits"

for i_sim,sim in enumerate(simulations):
    for i_ang, angle in enumerate(angles) :
        for i_im,imgr in enumerate(imagers):
            print(i_sim,i_ang,i_im)
            os.system(f"rm ../workspace/{i_sim}_{i_ang}_{i_im}_??.png")

            print(sim,imgr,angle)
            fname = file_base.format(skirt_data_dir=skirt_data_dir,sim=sim,imgr=imgr,angle=angle)
            with fits.open(fname) as f :
                d = f[0].data # image maps
                mic = f[1].data # wavelengths in microns

            d = np.log10(d)
            print(np.nanmax(d),np.nanmin(d[np.isfinite(d)]))
            nmic = d.shape[0]
            for imic in range(nmic):
                fig,sp = plt.subplots(figsize=(6,6),dpi=200)
                im = d[imic,:]
                sp.imshow(im,cmap='viridis')
                wl = mic[imic][0]
                sp.set_title(f"{wl:.2f}"+r"$\mu$m, angle = "+f"{angle}" + r"$^\circ$")
                fig.savefig(f"../workspace/{i_sim}_{i_ang}_{i_im}_{imic:02d}.png")
                plt.close(fig)
            cmd = f"ffmpeg -y -r 24 -i ../workspace/{i_sim}_{i_ang}_{i_im}_%02d.png -c:v mpeg4 -q:v 1 ../quick_out/{sim}_big_{imgr}_{angle}_total.mp4"
            # cmd = f"ffmpeg -y -r 24 -hide_banner -loglevel quiet -stats -i ../workspace/{i_sim}_{i_ang}_{i_im}_%02d.png -c:v mpeg4 -q:v 1 ../quickout/{sim}_big_{imgr}_{angle}_total.mp4"
            os.system(cmd)