import h5py
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import vapeplot
from scipy import interpolate
# vapeplot.set_palette("vaporwave")
plt.style.use('dark_background')

runs = [["norad_prod", ["small_ecc","small_circ"]],
        ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"]],
        ["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"]]
        ]

# norad_run = ["norad_prod", ["small_ecc","small_circ"]]
# rad_full =  ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"]]
# rad_earlier = ["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"]]

filebase = "../analysis_out/torque_{}_{}.hdf5"
filenames = [[filebase.format(x[0],x[1][0]),filebase.format(x[0],x[1][1])]
             for x in runs]

h5_files = [ [h5py.File(filename, "r") for filename in x] for x in filenames]
gizmoDir = "/srv/djw1g16/gizmos/"

# binary_key_bases = ["BH_acc_J","BH_acc_mom","BH_force","BH_grav_J","BH_grav_mom","BH_mass","BH_pos","BH_vel"]
# combined_keys = ["acc_torque","angmom","dist","grav_torque"]

binary_key_bases = ["BH_acc_J","BH_grav_J"]
combined_keys = ["acc_torque","grav_torque"]
t_offset = 0.0002128*0.9778e9/1.e6


ny = len(combined_keys) + len(binary_key_bases)
pix_width = 3000
inch_width = 6
inch_height = 4

cmap = cm.get_cmap("gnuplot")
colours = cmap(np.linspace(0.2,1,6))

def glow_plot(ax,x,y,c,n=5,label=None,**kwargs):
    ax.plot(x, y, c=c, lw=1, alpha=1., label=label, **kwargs)
    # if n>1:
    #     for width in np.linspace(2,10,n-1):
    #         ax.plot(x,y,c=c,lw=width,alpha=1./width,**kwargs)

t_offsets = [[0,None,0.0002128*0.9778e9],[0,None,0.0002128*0.9778e9]]

t_offsets[0][1] = np.array(h5_files[0][0]["time"])[-1]
t_offsets[1][1] = np.array(h5_files[0][1]["time"])[-1]

# t_offset = 0.0002128*0.9778e9/1.e6
t_orbit = 13.7e-3 # Myr
t_orbit *=1.e6 # yr

nx = 2
fig,sp = plt.subplots(ny,nx,figsize=(inch_width*nx, inch_height * ny), dpi=pix_width/inch_width, constrained_layout=True)
for isim in range(2):
    iy = 0
    ix = isim
    for key_basis in binary_key_bases:
        ax = sp[iy,ix]
        ax2 = sp[iy+2,ix]
        ax.set_ylabel(key_basis)
        for ibh,bh in enumerate(["_1","_2"]):
            key = key_basis+bh
            for irun,offset in enumerate(t_offsets[isim]):
                dt = t_offsets[isim][irun]
                time = np.array(h5_files[irun][isim]["time"]) + dt
                val = np.array(h5_files[irun][isim][key])
                if irun == 1 :
                    val += np.array(h5_files[0][isim][key])[-1]
                elif irun == 2 :
                    y = np.array(h5_files[0][isim][key])
                    x = np.array(h5_files[0][isim]["time"])
                    if len(y.shape) == 2 :
                        y = y[:, 2]
                    dy = interpolate.interp1d(x, y)(dt)
                    val += dy
                elif irun==0: # interpolate data to make timing consistent
                    big_dt = time[-1]-time[-2]
                    interp_t = np.arange(0, time[-1], big_dt)
                    if len(val.shape)==2:
                        y = val[:,2]
                    else:
                        y = val
                    val = interpolate.interp1d(time,y)(interp_t)
                    time = interp_t
                time/=t_orbit
                x = time
                if len(val.shape)==2:
                    y = val[:,2]
                else:
                    y = val
                glow_plot(ax, x, y, colours[ibh * 3+irun], label="z" + bh)

                dydx = np.gradient(y,x)
                glow_plot(ax2, x, dydx, colours[ibh * 3+irun], label="z" + bh)
        iy+=1
    sp[3,ix].set_yscale('symlog', linthresh=1e-7)
    sp[2,ix].set_yscale('symlog', linthresh=1e-8)
    # for key in combined_keys:
    #     ax = sp[iy,ix]
    #     ax.set_ylabel(key)
    #     for irun,offset in enumerate(t_offsets[isim]):
    #         dt = t_offsets[isim][irun]
    #         time = np.array(h5_files[irun][isim]["time"])+dt
    #         val = np.array(h5_files[irun][isim][key])
    #         if irun==1:
    #             val+=np.array(h5_files[0][isim][key])[-1]
    #         elif irun==2:
    #             y = np.array(h5_files[0][isim][key])
    #             x = np.array(h5_files[0][isim]["time"])
    #             if len(y.shape)==2:
    #                 y = y[:,2]
    #             dy = interpolate.interp1d(x,y)(dt)
    #             val+=dy
    #         elif irun == 0 :  # interpolate data to make timing consistent
    #             big_dt = time[-1] - time[-2]
    #             interp_t = np.arange(0, time[-1], big_dt)
    #             if len(val.shape) == 2 :
    #                 y = val[:, 2]
    #             else :
    #                 y = val
    #             val = interpolate.interp1d(time, y)(interp_t)
    #             time = interp_t
    #         time /= t_orbit
    #         if len(val.shape)==2:
    #             glow_plot(ax, time, val[:, 2], colours[irun], label="z" + bh)
    #         else:
    #             glow_plot(ax, time, val, colours[irun], label="z" + bh)
    #     ax.set_yscale('symlog',linthresh=1.e-13)
    #     iy+=1

fig.savefig(f"../quick_out/all_torque_cont.pdf")

# nx = len(runs)
# inch_size = 4
# fig, sp = plt.subplots(ny, nx, figsize=(inch_size*nx, inch_size * ny), dpi=pix_width/inch_size, constrained_layout=True)#, sharey='row')
#
#
# for ix,(run_id,output_dir) in enumerate(runs):
#     print(run_id,output_dir)
#     filename = f"../analysis_out/torque_{run_id}_{output_dir}.hdf5"
#     f = h5py.File(filename, "r")
#     iy = 0
#     time = f["time"]
#
#     sp[0,ix].set_title(f"{run_id}_{output_dir}")
#
#     for key_basis in binary_key_bases:
#         ax = sp[iy,ix]
#         ax.set_ylabel(key_basis)
#         for ibh,bh in enumerate(["_1","_2"]):
#             key = key_basis+bh
#             val = np.array(f[key])
#             if len(val.shape)==2:
#                 # vector, plot all 3
#                 glow_plot(ax, time, val[:, 2], colours[ibh * 3], label="z" + bh)
#                 # for idir in range(3):
#                 #     glow_plot(ax,time,val[:,idir],colours[idir+ibh*3],label="xyz"[idir]+bh)
#                     # ax.plot(time,val[:,idir],label="xyz"[idir]+bh)
#             else:
#                 glow_plot(ax, time, val, colours[ibh * 3], label=bh)
#                 # ax.plot(time,val,label=bh)
#         ax.legend()
#         iy+=1
#
#     for key in combined_keys:
#         ax = sp[iy,ix]
#         ax.set_ylabel(key)
#         val = np.array(f[key])
#         if len(val.shape) == 2 :
#             # vector, plot all 3
#             glow_plot(ax, time, val[:, 2], colours[0], label="z")
#             # for idir in range(3) :
#             #     # ax.plot(time, val[:, idir], label="xyz"[idir])
#             #     glow_plot(ax, time, val[:, idir], colours[idir], label="xyz"[idir])
#             #     ax.legend()
#         else :
#             glow_plot(ax, time, val, colours[0])
#             # ax.plot(time, val)
#         iy+=1
#
#     sp[3,ix].set_ylim(-2e-12,2e-12)
#     sp[2,ix].set_ylim(-2e-12,2e-12)
#     f.close()
#
# fig.savefig(f"../quick_out/all_torque.pdf")
# fig.savefig(f"../quick_out/all_torque.png")
# df = pd.read_hdf(f"../analysis_out/torque_{run_id}_{output_dir}.hdf5")
    # print(df)