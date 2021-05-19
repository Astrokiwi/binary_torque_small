# import scipy
# from scipy import signal
import numpy as np
import h5py
from scipy import optimize
from functools import partial

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

def constant_to_straight(x,a,m,xc=0.1):
    return np.piecewise(
        x,
        [x<xc, x>=xc],
        [ lambda x: a, lambda x: a-m*(x-xc)]
        )

def constant(x,a):
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

run_labels = ["ecc","circ"]
set_labels = ["norad","rad","rad_earlier"]

fout = h5py.File("/srv/djw1g16/paper_devel/binary_torque_small/analysis_out/power.hdf5",'w')

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

                # ax = sp[iset, irun + 2*ikey]
                pos_slice = (x>0)
                pos_paua = paua[idx][pos_slice]
                x = x[pos_slice]

                logx = np.log10(x)
                logy = np.log10(pos_paua)

                popt,pcov = find_breakpoint(split_poly,logx,logy)
                perr = np.sqrt(np.diag(pcov))

                iy = ikey * 3 + iset
                ix = irun

                fout[f"power_{ix}_{iy}_{ibh}"] = pos_paua
                fout[f"freq_{ix}_{iy}_{ibh}"] = x
                fout[f"popt_{ix}_{iy}_{ibh}"] = popt
                fout[f"perr_{ix}_{iy}_{ibh}"] = perr


fout.close()
f.close()
