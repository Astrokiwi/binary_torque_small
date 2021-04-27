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

def split_powerlaw(x,a,m,m2,xc=0.1):
    return np.piecewise(
        x,
        [x<xc, x>=xc],
        [ lambda x: a*(x/xc)**m2, lambda x: a*(x/xc)**m]
        )

def split_poly(x,a,m,m2,xc=0.1):
    return np.piecewise(
        x,
        [x<xc, x>=xc],
        [ lambda x: a+(x-xc)*m2, lambda x: a+(x-xc)*m]
        )


def constant_to_powerlaw(x,a,m,xc=0.1):
    return np.piecewise(
        x,
        [x<xc, x>=xc],
        [ lambda x: a, lambda x: a*(x/xc)**m]
        )
    # y1 = (a*x**0)[x<xc]
    # y2 = (a*(x/xc)**c)[x>=xc]
    # y = np.concatenate((y1,y2))
    # return y
    # return a/(1+b*x**(-c))

def constant_to_straight(x,a,m,xc=0.1):
    return np.piecewise(
        x,
        [x<xc, x>=xc],
        [ lambda x: a, lambda x: a-m*(x-xc)]
        )
    # y1 = (a*x**0)[x<xc]
    # y2 = (a-m*(x-xc))[x>=xc]
    # y = np.concatenate((y1,y2))
    # return y

def constant(x,a):
    return a*np.ones(x.shape)

def powerlaw(x,m):
    return a*np.ones(x.shape)


def find_breakpoint(f,x,y,n=100):
    prams = None
    best_fit = 1.e20
    best_pcov = None
    # xc_range = np.logspace(np.log10(x[0]),np.log10(x[-1]),n)
    xc_range = np.linspace(x[0],x[-1],n)
    for xc in xc_range:
        f_xc = partial(f,xc=xc)
        popt, pcov = optimize.curve_fit(f_xc,x,y,p0=[y[0],-2,0],bounds=([-np.inf,-5,-5],[np.inf,1,1]))
        perr = np.sum(np.abs(y - f_xc(x, *popt)))
        if perr<best_fit:
            prams = np.append(popt,xc)
            best_fit = perr
            best_pcov = pcov
    return prams,best_pcov


f=h5py.File("/srv/djw1g16/paper_devel/binary_torque_small/analysis_out/tidy_torque.hdf5",'r')

fig, sp = plt.subplots(3, 4, constrained_layout=True, figsize=(12, 6), dpi=200,sharex=True)
for irun in range(2) :
    for iset in range(3) :
        for ibh in range(2) :
            for ikey,key in enumerate(["acc","grav"]):
                y = np.array(f[f"BH_{key}_J_{ibh + 1}_{iset}_{irun}"])
                t = np.array(f[f"time_{iset}_{irun}"])
                # cut last one, where dt is not consistent
                t = t[:-1]
                y = y[:-1]

                dydt = np.gradient(y, t)

                dt = t[-1] - t[-2]

                paua = np.abs(np.fft.fft(dydt)) ** 2
                freqs = np.fft.fftfreq(dydt.size, dt)
                idx = np.argsort(freqs)
                x = freqs[idx]

                ax = sp[iset, irun + 2*ikey]
                pos_slice = (x>0)
                pos_paua = paua[idx][pos_slice]
                x = x[pos_slice]

                ax.loglog(x, pos_paua)

                logx = np.log10(x)
                logy = np.log10(pos_paua)
                # popt,pcov = optimize.curve_fit(constant_to_straight,logx,logy)
                # print(popt)
                # ax.loglog(10**logx,10**constant_to_straight(logx,*popt))

                # popt,pcov = optimize.curve_fit(constant_to_powerlaw,x,pos_paua,p0=[pos_paua[0],0.1,-2])
                # try:
                popt,pcov = find_breakpoint(split_poly,logx,logy)
                perr = np.sqrt(np.diag(pcov))
                # print(pcov,perr)
                # print(popt)
                # print(pcov)
                # ax.loglog(x,split_powerlaw(x,*popt))
                label_txt = rf"$\alpha={popt[2]:.1f}\pm{perr[2]:.1f},{popt[1]:.1f}\pm{perr[1]:.1f}$"
                # print(label_txt)
                ax.loglog(10**logx,10**split_poly(logx,*popt),ls='--',lw=1,zorder=2.1,label=label_txt)
                # except RuntimeError as e:
                #     print(str(e))
                #     print("Runtime Error - carrying on")
                # ax.set_ylim(pos_paua[0]*0.9,pos_paua[-1]*1.1)
                # if ibh==0:
                #     ax.loglog(x,x**-1)
                #     ax.loglog(x,x**-2)
                #     ax.loglog(x,np.ones(x.size)*pos_paua[0])

for ix in range(4):
    for iy in range(3):
        ax = sp[iy,ix]
        ax.set_xticks([0.01,0.05,0.1,0.5,1.,5.,10.])
        ax.legend(fontsize='xx-small')

for ix in range(4):
    sp[-1,ix].set_xlabel(r'$\omega$ ($1/T$)')
for iy in range(3):
    sp[iy,0].set_ylabel(r'$P$')
fig.savefig("../figures_out/powerspectra.pdf")

f.close()
