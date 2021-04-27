import h5py
import numpy as np
import pandas as pd
from scipy import interpolate

runs = [["norad_prod", ["small_ecc","small_circ"]],
        ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"]],
        ["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"]]
        ]

filebase = "../analysis_out/torque_{}_{}.hdf5"
filenames = [[filebase.format(x[0],x[1][0]),filebase.format(x[0],x[1][1])]
             for x in runs]

h5_files = [ [h5py.File(filename, "r") for filename in x] for x in filenames]
gizmoDir = "/srv/djw1g16/gizmos/"

binary_key_bases = ["BH_acc_J","BH_grav_J"]
combined_keys = ["acc_torque","grav_torque"]
t_offset = 0.0002128*0.9778e9/1.e6

t_offsets = [[0,None,0.0002128*0.9778e9],[0,None,0.0002128*0.9778e9]]

t_offsets[0][1] = np.array(h5_files[0][0]["time"])[-1]
t_offsets[1][1] = np.array(h5_files[0][1]["time"])[-1]

# t_offset = 0.0002128*0.9778e9/1.e6
t_orbit = 13.7e-3 # Myr
t_orbit *=1.e6 # yr

# calculate angular momentum for normalisation
semi_major = 0.035 # pc
eccentricities = [0.5,0.] # unitless
masses = [[500,1500],[1000,1000]] # Msun
G = 4.5e-15 # pc^3 /year**2 /Msun

angular_momentum = [
    mass[0]*mass[1]*np.sqrt(G*semi_major*(1-eccentricity**2)/(mass[0]+mass[1]))

    for eccentricity,mass in zip(eccentricities,masses)
    ]

f_out = h5py.File("../analysis_out/tidy_torque.hdf5", "w")

for isim in range(2):
    for key_basis in binary_key_bases:
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
                    big_dt = time[-2]-time[-3]
                    interp_t = np.arange(0, time[-1], big_dt)
                    if len(val.shape)==2:
                        y = val[:,2]
                    else:
                        y = val
                    val = interpolate.interp1d(time,y)(interp_t)
                    time = interp_t
                time/=t_orbit
                print(time)
                x = time
                if len(val.shape)==2:
                    y = val[:,2]
                else:
                    y = val
                h5_key = f"/{key}_{irun:01d}_{isim:01d}"
                # normalize by J_initial
                f_out[h5_key] = y/angular_momentum[isim]
                h5_time_key = f"/time_{irun:01d}_{isim:01d}"
                if h5_time_key not in f_out:
                    f_out[h5_time_key] = x

f_out.close