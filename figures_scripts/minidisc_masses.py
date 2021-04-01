import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from scipy import interpolate

data_dir = "../analysis_out/"

norad_folder = "norad_prod"
rad_folder = "rad_prod"
sim_name_bases = ["small_ecc","small_circ"]
rad_format = "rad_{}_full"
rad_format_earlier = "rad_{}_earlier"

data_file_format = data_dir+"disc_components_{}_{}.txt"
t_orbit = 13.7e-3 # Myr
m_bh = 2.e3

t_offset = 0.0002128*0.9778e9/1.e6
norad_data_unprocessed = [pd.read_csv(data_file_format.format(norad_folder,sim_name_bases[irun]),sep=' ') for irun in [0,1]]
rad_data = [pd.read_csv(data_file_format.format(rad_folder,rad_format.format(sim_name_bases[irun])),sep=' ') for irun in [0,1]]
rad_data_earlier = [pd.read_csv(data_file_format.format(rad_folder,rad_format_earlier.format(sim_name_bases[irun])),sep=' ') for irun in [0,1]]

# interpolate norad_data to make timing consistent
norad_data = []
for df in norad_data_unprocessed:
    dt = df["t"].iloc[-1]-df["t"].iloc[-2]
    interp_t = np.arange(0,df["t"].iloc[-1],dt)

    new_df = pd.DataFrame()
    new_df["t"] = interp_t
    new_df["m1"] = interpolate.interp1d(df["t"],df["m1"])(interp_t)
    new_df["m2"] = interpolate.interp1d(df["t"],df["m2"])(interp_t)
    norad_data.append(new_df)

# yticks = [0,5e-4,1e-3,3e-3,5e-3,1e-2]

fig,sp = plt.subplots(1,2,figsize=(9,3),sharex=True,sharey=True,constrained_layout=True)
for irun in range(2):
    for idisc in [1,2]:
        for t,m,label in [(norad_data[irun]["t"], norad_data[irun][f"m{idisc}"], f"disc{idisc}"),
                          (rad_data[irun]["t"] + norad_data[irun]["t"].iloc[-1], rad_data[irun][f"m{idisc}"], None),
                          (rad_data_earlier[irun]["t"] + t_offset, rad_data_earlier[irun][f"m{idisc}"], None)
                          ]:
            sp[irun].plot(t / t_orbit, m / m_bh, label=label)
            # orbit_t = t/t_orbit
            # interp_t = np.arange(int(orbit_t[0]),orbit_t[-1],1)
        # sp[irun].plot(norad_data[irun]["t"]/t_orbit,norad_data[irun][f"m{idisc}"]/m_bh,label=f"disc{idisc}")
        # sp[irun].plot((rad_data[irun]["t"]+norad_data[irun]["t"].iloc[-1])/t_orbit,rad_data[irun][f"m{idisc}"]/m_bh)
        # sp[irun].plot((rad_data_earlier[irun]["t"]+t_offset)/t_orbit,rad_data_earlier[irun][f"m{idisc}"]/m_bh)
    # sp[irun].legend()
    sp[irun].minorticks_on()
    sp[irun].xaxis.set_minor_locator(plt.MultipleLocator(1))
    sp[irun].set_xlim([0,None])
    # sp[irun].set_ylim([0,1.e-2])
    sp[irun].set_ylim([0,0.003])
    # sp[irun].set_yscale('symlog',linthresh=1.e-3)
    # sp[irun].set_yticks(yticks)
    # sp[irun].set_yticklabels(yticks)
    sp[irun].set_xlabel(r"$t$ (orbits)")
sp[0].set_ylabel(r"$M_d/(M_1+M_2)$")
sp[0].set_title("Elliptical")
sp[1].set_title("Circular")
plt.savefig("../figures_out/inner_mass.pdf")