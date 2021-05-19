# import scipy
# from scipy import signal
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
from scipy import optimize
from scipy import interpolate
from functools import partial
import sys

# def constant_to_powerlaw(x,a,b,c):
#     return a/(1+b*x**(-c))

def split_poly(x,a,m,m2,xc=0.1):
    return np.piecewise(
        x,
        [x<xc, x>=xc],
        [ lambda x: a+(x-xc)*m2, lambda x: a+(x-xc)*m]
        )

def format_with_uncertainty(x,e,precision=2):
    round_factor = -np.floor(np.log10(e)).astype(np.int) + precision-1
    x_round = np.round(x,round_factor)
    e_round = np.round(e,round_factor)

    if round_factor<0:
        strout = fr"{{{x_round:.0f}}} \pm {{{e_round:.0f}}}"
    else:
        strout_0 = fr"{{0:.{round_factor:0d}f}} \pm {{1:.{round_factor:0d}f}}"
        strout = strout_0.format(x_round,e_round)

    strout = r"$"+strout+r"$"
    return strout

f=h5py.File("/srv/djw1g16/paper_devel/binary_torque_small/analysis_out/power.hdf5",'r')

run_labels = ["ecc","circ"]
set_labels = ["norad","rad","rad_earlier"]

nx = 2
ny = 6
fig, sp = plt.subplots(ny+1, nx, constrained_layout=True, figsize=(nx*3, ny*2.), dpi=200,sharex=True, gridspec_kw={"height_ratios":[1,1,1,0.2,1,1,1]})

y_plots = [0,1,2,4,5,6]

for irun in range(2) :
    for iset in range(3) :
        for ikey,key in enumerate(["acc","grav"]):
            # ax = sp[iset, irun + 2 * ikey]
            ax = sp[ikey*3+iset, irun]
            ax.set_title(f"{run_labels[irun]}_{set_labels[iset]}_{key}")

colour_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for irun in range(2) :
    for iset in range(3) :
        for ikey,key in enumerate(["acc","grav"]):
            y_mins = []
            y_maxes = []

            iy = ikey * 4 + iset
            ix = irun
            ax = sp[iy, ix]

            for ibh in range(2) :
                iy = ikey * 3 + iset

                pos_paua = f[f"power_{ix}_{iy}_{ibh}"]
                x = f[f"freq_{ix}_{iy}_{ibh}"]
                popt = f[f"popt_{ix}_{iy}_{ibh}"]
                perr = f[f"perr_{ix}_{iy}_{ibh}"]

                ax.loglog(x, pos_paua, c=colour_cycle[iset+ibh*3])

                logx = np.log10(x)
                logy = np.log10(pos_paua)

                label_txt = rf"$\alpha={popt[2]:.1f}\pm{perr[2]:.1f},{popt[1]:.1f}\pm{perr[1]:.1f}$"
                label_txt = format_with_uncertainty(popt[2],perr[2])+","+format_with_uncertainty(popt[1],perr[1])
                # ax.loglog(10**logx,10**split_poly(logx,*popt),ls='--',lw=1,zorder=2.1,label=label_txt,c=colour_cycle[7+ibh])
                ax.loglog(10**logx,10**split_poly(logx,*popt),label=label_txt,c=colour_cycle[7+ibh])

                y_mins.append(np.min(pos_paua))
                y_maxes.append(np.max(pos_paua))

            y0 = np.min(y_mins)
            y1 = np.max(y_maxes)
            dy = y1/y0
            ax.set_ylim(y1/(dy**1.4),y0*(dy**1.1))

for ix in range(nx) :
    sp[3,ix].remove()
    for iy in y_plots :
        ax = sp[iy,ix]
        ticks = np.array([0.01,0.05,0.1,0.5,1.,5.,10.])
        # labels = 1./ticks
        ax.set_xticks(ticks)
        # ax.set_xticklabels(["100","20","10","2","1","0.2","0.1"])

        ax.legend(fontsize='small')

for ix in range(nx):
    sp[-1,ix].set_xlabel(r'$f$ ($1/T$)')
    # sp[-1,ix].set_xlabel(r'Period ($T$)')
for iy in y_plots:
    sp[iy,0].set_ylabel(r'$P$')
fig.savefig("../figures_out/powerspectra.pdf")

f.close()
