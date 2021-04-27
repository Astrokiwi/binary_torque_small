import numpy as np
import h5py
from num2tex import num2tex

f = h5py.File("../analysis_out/tidy_torque.hdf5",'r')

def get_exponent(x):
    return np.floor(np.log10(np.abs(x))).astype(np.int)

def format_with_uncertainty(x,e,precision=2):
    fmt = fr"$({{0:.{precision}f}} \pm {{1:.{precision}f}}) \times 10^{{{{{{2}}}}}}$"

    exponents = np.minimum(get_exponent(x),get_exponent(e))
    bases = x/10.**exponents
    error_bases = e/10.**exponents
    # TODO: implement for arrays?
    return fmt.format(bases,error_bases,exponents)


# UNITS ??
for plot_key in ["grav","acc"]:
    mu = np.zeros((2,2,3))
    sd = np.zeros((2,2,3))
    for isim in range(2):
        for bh in [1,2] :
            for irun in range(3):
                time = np.array(f[f"time_{irun}_{isim}"])

                h5_key = f"BH_{plot_key}_J_{bh}_{irun}_{isim}"

                angmom = -np.array(f[h5_key])
                torque = np.gradient(angmom,time)

                mu[isim,bh-1,irun] = np.mean(torque)
                sd[isim,bh-1,irun] = np.std(torque)/np.sqrt(len(torque))
    for irun in range(3):
        line = []
        for isim in range(2) :
            for bh in range(2) :
                line+=[format_with_uncertainty(mu[isim,bh,irun],sd[isim,bh,irun])]
                # line+=[fr"${num2tex(mu[isim,bh,irun],precision=2):.2g} \pm {num2tex(sd[isim,bh,irun],precision=2):.0g}$"]
        print(r" & ".join(line)+r"\\")
    print(r"\\")
