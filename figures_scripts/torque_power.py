# import scipy
# from scipy import signal
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
from scipy import optimize
from scipy import signal
from scipy import interpolate
from functools import partial
import sys

# def constant_to_powerlaw(x,a,b,c):
#     return a/(1+b*x**(-c))

def exp_drop(x,a,b):
    return a*np.exp(-b*x)#+1.e-10

# def exp_drop(x,a,m,m2,xc):
#     return np.piecewise(
#         x,
#         [x<xc, x>=xc],
#         [ lambda x: a+(x-xc)*m2, lambda x: a+(x-xc)*m]
#         )


# def exp_drop(x,m,c):
#     return -m*x+c


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

f=h5py.File("../analysis_out/power.hdf5",'r')

run_labels = ["ecc","circ"]
set_labels = ["norad","rad","rad_earlier"]

nx = 2
ny = 6
scale = 1.5
# fig, sp = plt.subplots(ny+1, nx, constrained_layout=True, figsize=(nx*3*scale, ny*2.*scale), dpi=200,sharex='col', gridspec_kw={"height_ratios":[1,1,1,0.2,1,1,1]})
fig, sp = plt.subplots(2+1, nx, constrained_layout=True, figsize=(nx*3*scale, 2*2.*scale), dpi=200,sharex='col', gridspec_kw={"height_ratios":[1,0.2,1]})

fig_torque, sp_torque = plt.subplots(ny+1, nx, constrained_layout=True, figsize=(nx*3.*scale, ny*2.*scale), dpi=200,gridspec_kw={"height_ratios":[1,1,1,0.2,1,1,1]})

y_plots = [0,1,2,4,5,6]

for irun in range(2) :
    for iset in range(3) :
        for ikey,key in enumerate(["acc","grav"]):
            # ax = sp[iset, irun + 2 * ikey]
            iy = ikey * 4 + iset
            ax_abs = sp_torque[iy, irun]
            ax_abs.set_title(f"{run_labels[irun]}_{set_labels[iset]}_{key}")


for irun in range(2) :
    for ikey, key in enumerate(["acc", "grav"]) :
        iy = ikey * 2
        ax_pow = sp[iy, irun]
        ax_pow.set_title(f"{run_labels[irun]}_{set_labels[iset]}_{key}")

colour_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for irun in range(2) :
    for ikey,key in enumerate(["acc","grav"]):
        y_mins = []
        y_maxes = []

        # iy = ikey * 4 + iset
        iy = ikey * 2
        ix = irun
        ax_pow = sp[iy, ix]

        for ibh in range(2) :
            # iy = ikey * 3 + iset
            iy = ikey * 3

            pos_paua = f[f"power_{ix}_{iy}_{ibh}"]
            x = f[f"freq_{ix}_{iy}_{ibh}"]
            popt = f[f"popt_{ix}_{iy}_{ibh}"]
            perr = f[f"perr_{ix}_{iy}_{ibh}"]

            ax_pow.loglog(x, pos_paua, c=colour_cycle[ibh*3])

            logx = np.log10(x)
            logy = np.log10(pos_paua)

            label_txt = rf"$\alpha={popt[2]:.1f}\pm{perr[2]:.1f},{popt[1]:.1f}\pm{perr[1]:.1f}$"
            label_txt = format_with_uncertainty(popt[2],perr[2])+","+format_with_uncertainty(popt[1],perr[1])
            # ax.loglog(10**logx,10**split_poly(logx,*popt),ls='--',lw=1,zorder=2.1,label=label_txt,c=colour_cycle[7+ibh])
            ax_pow.loglog(10**logx,10**split_poly(logx,*popt),label=label_txt,c=colour_cycle[7+ibh])

            y_mins.append(np.min(pos_paua))
            y_maxes.append(np.max(pos_paua))


        y0 = np.min(y_mins)
        y1 = np.max(y_maxes)
        dy = y1/y0
        ax_pow.set_ylim(y1/(dy**1.4),y0*(dy**1.1))

for irun in range(2) :
    for iset in range(3) :
        for ikey,key in enumerate(["acc","grav"]):
            y_mins = []
            y_maxes = []

            iy = ikey * 4 + iset
            ix = irun
            ax_abs = sp_torque[iy, ix]

            for ibh in range(2) :
                iy = ikey * 3 + iset
                torque = f[f"tau_{ix}_{iy}_{ibh}"]
                time = np.array(f[f"t_{ix}_{iy}_{ibh}"])

                t = time-time[0]
                ax_abs.plot(t,torque, c=colour_cycle[iset+ibh*3])

                # if iset==0:
                #     fit_slice = (t>20.)
                #     t_fit = t[fit_slice]
                #     torque_fit = np.abs(torque[fit_slice])
                #     popt, pcov = optimize.curve_fit(exp_drop,t_fit,torque_fit)
                #     # print(popt)
                #     # ax_abs.plot(t,exp_drop(t,*popt), c=colour_cycle[iset+ibh*3],ls='--')
                #     ax_abs.plot(t[fit_slice],torque[fit_slice]/exp_drop(t[fit_slice],*popt), c=colour_cycle[iset+ibh*3])


            if iset == 0 and ikey==0:
                if irun==0:
                    ax_abs.set_yscale('symlog',linthresh=1e-5)
                elif irun==1:
                    ax_abs.set_yscale('symlog', linthresh=1e-5)


for ix in range(nx) :
    # sp[3,ix].remove()
    sp[1,ix].remove()
    sp_torque[3,ix].remove()

for ix in range(nx) :
    # for iy in y_plots :
    for iy in [0,2] :
        ax = sp[iy,ix]
        ticks = np.array([0.01,0.05,0.1,0.5,1.,5.,10.])
        # labels = 1./ticks
        ax.set_xticks(ticks)
        # ax.set_xticklabels(["100","20","10","2","1","0.2","0.1"])

        ax.legend(fontsize='small')

for ix in range(nx):
    sp[-1,ix].set_xlabel(r'$f$ ($1/T$)')
    sp_torque[-1,ix].set_xlabel(r'$t$ ($T$)')
    # sp[-1,ix].set_xlabel(r'Period ($T$)')
for iy in [0,2]:
    sp[iy,0].set_ylabel(r'$P$')
for iy in y_plots :
    sp_torque[iy, 0].set_ylabel(r'$\tau$')
fig.savefig("../figures_out/powerspectra.pdf")
fig_torque.savefig("../figures_out/torque_time.pdf")

f.close()
