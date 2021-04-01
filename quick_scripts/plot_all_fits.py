import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from astropy.io import fits
import numpy as np

simulations = ["binary","unary"]
# imagers = ["K","L","12"]
imagers = ["K","bands"]
# angles = ["face","nearface","edge","nearedge"]
angles = [f"{inc:02d}" for inc in np.linspace(0,90,11).astype(np.int)]

nx = len(imagers)*len(simulations)
ny = len(angles)

with open("../local_config/skirt_data_dir") as f:
    skirt_data_dir = f.readline()

# file_base = "{skirt_data_dir}/{sim}_lowres_{imgr}_{angle}_total.fits"
# file_base = "{skirt_data_dir}/{sim}_zoom_{imgr}_{angle}_total.fits"
file_base = "{skirt_data_dir}/big_data/{sim}_big_{imgr}_{angle}_total.fits"

# plot flattened
fig,sp = plt.subplots(ny,nx,figsize=(nx*4,ny*3),dpi=100,constrained_layout=True,sharex=True,sharey=True)
for i_sim,sim in enumerate(simulations):
    for i_im,imgr in enumerate(imagers):
        for i_ang,angle in enumerate(angles):
            print(sim,imgr,angle)
            fname = file_base.format(skirt_data_dir=skirt_data_dir,sim=sim,imgr=imgr,angle=angle)
            with fits.open(fname) as f :
                d = f[0].data
            d = np.sum(d, axis=0) # flatten across wavelengths
            d = np.log10(d)

            ix = i_im + i_sim * len(imagers)
            iy = i_ang
            ax = sp[iy,ix]
            # ax.imshow(d,norm=colors.LogNorm())
            ax.imshow(d)

print("Saving")
fig.savefig("../quick_out/flattened_fits_big_all.png")
print("Saved")

plt.close(fig)