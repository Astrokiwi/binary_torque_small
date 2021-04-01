# Import libraries to do our magic

import os
this_dir, this_filename = os.path.split(__file__)


import numpy as np

from gtools import gizmo_tools
from multiprocessing import Pool

from scipy import linalg
import functools

import matplotlib as mpl

mpl.use('Agg')

from scipy import interpolate

import tqdm
import deepdish as dd

runs = [["norad_prod", "small_ecc"],
        ["norad_prod", "small_circ"],
        ["rad_prod","rad_small_ecc_full"],
        ["rad_prod","rad_small_circ_full"],
        ["rad_prod","rad_small_ecc_earlier"],
        ["rad_prod","rad_small_circ_earlier"],
        ]

gizmoDir = "/srv/djw1g16/gizmos/"


def load_gadget(run_id, output_dir, snap_str):
    return gizmo_tools.load_gizmo_nbody(run_id, output_dir, snap_str=snap_str, gizmoDir=gizmoDir
                                        , load_binary_headers=True, only_header=True)


def extract_bh_data(run_id, output_dir, gizmodir=gizmoDir, snap0=0, maxsnapf=-1, nprocs=64):
    snapi = 0
    if snap0 > 0:
        snapi = snap0
    snapf = gizmo_tools.lastConsecutiveSnapshot(run_id, output_dir, False, gizmoDir=gizmodir)

    # max run
    if -1 < maxsnapf < snapf:
        print("Forcing snapf from {} to {}".format(snapf, maxsnapf))
        snapf = maxsnapf

    print("nfiles:", snapf - snapi + 1)

    isnaps = range(snapi, snapf + 1)

    snap_strs = [f"{i:03d}" for i in isnaps]
    extract_header_func = functools.partial(load_gadget, run_id, output_dir)

    with Pool(processes=nprocs) as pool:
        bh_binary_alltimes = pool.map(extract_header_func, snap_strs)

    # convert list of dicts of 1D arrays to dict of 2D arrays
    d = {k: np.array([dic[k] for dic in bh_binary_alltimes]) for k in bh_binary_alltimes[0]}

    return d


def calc_torque_etc(d,resample=False):
    d['time']*=1.e6 # from Myr to years

    # fix BH vels
    dtime = np.gradient(d['time'])
    d['dist'] = linalg.norm(d['BH_pos_1'] - d['BH_pos_2'], axis=1)

    if resample:
        dt = np.max(dtime)
        tmax = np.max(d['time'])
        print(f"Resampling to dt {dt} years")
        d_new = dict()
        d_new['time'] = np.arange(0,tmax,dt)
        for key in d:
            d_new[key] = interpolate.interp1d(d["time"],d[key],fill_value="extrapolate")(d_new['time'])
        d = d_new


    d['BH_vel_1'] = np.gradient(d['BH_pos_1'], axis=0) / dtime[:, np.newaxis] / 1.e6
    d['BH_vel_2'] = np.gradient(d['BH_pos_2'], axis=0) / dtime[:, np.newaxis] / 1.e6

    d['angmom'] = d['BH_mass_1'][:, np.newaxis] * np.cross(d['BH_pos_1'], d['BH_vel_1']) \
                  + d['BH_mass_2'][:, np.newaxis] * np.cross(d['BH_pos_2'], d['BH_vel_2'])

    d['acc_torque'] = np.gradient(d['BH_acc_J_1']+d['BH_acc_J_2'],axis=0)/dtime[:,np.newaxis]
    d['grav_torque'] = np.gradient(d['BH_grav_J_1']+d['BH_grav_J_2'],axis=0)/dtime[:,np.newaxis]


def extract_save_torque(run_id, output_dir, **kwargs):
    # print(run_id,output_dir)
    d = extract_bh_data(run_id, output_dir, **kwargs)
    calc_torque_etc(d)
    dd.io.save(f"../analysis_out/torque_{run_id}_{output_dir}.hdf5",d)

if __name__ == '__main__':
    for run in tqdm.tqdm(runs):
        extract_save_torque(*run)
