import numpy as np
import h5py
from num2tex import num2tex

f = h5py.File("../analysis_out/tidy_torque.hdf5",'r')

def get_exponent(x):
    return np.floor(np.log10(np.abs(x))).astype(np.int)

def format_with_uncertainty(x,e,precision=1,significance=3):
    # fmt = fr"({{0:.{precision}f}} \pm {{1:.{precision}f}}) \times 10^{{{{{{2}}}}}}"
    #
    # exponents = np.minimum(get_exponent(x),get_exponent(e))
    # bases = x/10.**exponents
    # error_bases = e/10.**exponents
    # strout = fmt.format(bases,error_bases,exponents)

    round_factor = -np.floor(np.log10(e)).astype(np.int) + precision-1
    x_round = np.round(x,round_factor)
    e_round = np.round(e,round_factor)

    if round_factor<0:
        strout = fr"{{{x_round:.0f}}} \pm {{{e_round:.0f}}}"
    else:
        strout_0 = fr"{{0:.{round_factor:0d}f}} \pm {{1:.{round_factor:0d}f}}"
        strout = strout_0.format(x_round,e_round)
    # strout = fr"{{{x:.3f}}} \pm {{{e:.3f}}}"


    # if np.abs(bases)>np.abs(error_bases):
    if np.abs(x) > significance*np.abs(e) :
            strout = r"$\mathbf{"+strout+r"}$"
    else:
        strout = r"$"+strout+r"$"
    # TODO: implement for arrays?
    return strout


stage_name = [r"norad\_all", r"norad\_eqm", r"rad\_early", r"rad\_late"]

torque_name = {"grav":r"$\tau_g$",
                "acc":r"$\tau_a$"
    }

run_indices = [0,0,1,2]
run_offsets = [0,20,0,0]
nruns = len(run_indices)
assert nruns==len(run_offsets)

# UNITS = L0/T (set in tidy_up_torque)
for plot_key in ["acc","grav"]:
    mu = np.zeros((2,2,nruns))
    sd = np.zeros((2,2,nruns))
    for isim in range(2):
        for bh in [1,2] :
            for irun,run_index in enumerate(run_indices):
                time = np.array(f[f"time_{run_index}_{isim}"])

                h5_key = f"BH_{plot_key}_J_{bh}_{run_index}_{isim}"

                angmom = -np.array(f[h5_key])
                torque = np.gradient(angmom,time)

                run_offset = run_offsets[irun]
                if run_offset>0:
                    torque = torque[(time>run_offset)]

                mu[isim,bh-1,irun] = np.mean(torque)
                sd[isim,bh-1,irun] = np.std(torque)/np.sqrt(len(torque))


    for irun in range(nruns):
        line = []
        for isim in range(2) :
            for bh in range(2) :
                line+=[format_with_uncertainty(mu[isim,bh,irun]*1e6,sd[isim,bh,irun]*1e6,precision=2)]
                # line+=[fr"${num2tex(mu[isim,bh,irun],precision=2):.2g} \pm {num2tex(sd[isim,bh,irun],precision=2):.0g}$"]
        if irun==1:
            key_str = torque_name[plot_key]+" & "
        else:
            key_str = "~ & "
        print(key_str+stage_name[irun]+" & "+r" & ".join(line)+r"\\")
    print(r"\hline")
