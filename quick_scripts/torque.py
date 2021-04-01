# Import libraries to do our magic

import os
this_dir, this_filename = os.path.split(__file__)


import numpy as np

import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt

def next_axis(sp):
    irow = 0
    icol = 0

    while True:
        print(irow, icol)
        yield sp[irow, icol]
        icol += 1
        if icol >= sp.shape[1]:
            icol = 0
            irow += 1


def plot_torque(run_id, output_dir, d, tmax=0.05):
    time = d['time']
    dist = d['dist']
    angmom = d['angmom']

    plots = [
        #                'acc_force_1'
        #                 ,'grav_force_1'
        #                 ,'acc_force_2'
        #                 ,'grav_force_2'
        # 'bh_pos_1'
        # ,'bh_pos_2'
        # ,'gtorque'
        # , 'acc_torque'
        # ,'gcumtorque'
        # , 'acc_cumtorque'
        "acc_J"
        , "grav_J"
        , "acc_torque"
        , "grav_torque"
    ]

    nplots = len(plots)
    nrows = int(np.ceil(nplots / 2))
    ncols = 2

    scale = 3.
    fig, sp = plt.subplots(nrows, ncols, sharex=True, constrained_layout=True,
                           figsize=(ncols * scale * 2, nrows * scale), squeeze=False)

    print(plots, nplots,nrows,ncols,sp.shape)

    for iaxis, axis in enumerate(["x", "y", "z"]):
        na = next_axis(sp)

        if 'dist' in plots:
            ax = next(na)
            if iaxis == 0:
                ax.plot(time, dist)
                ax.set_ylabel(r'$d$ (pc)')

        if 'angmom' in plots:
            ax = next(na)
            ax.plot(time, angmom[:, iaxis], label=r"$L_{{{}}}$".format(axis))
            ax.set_ylabel('Angular momentum\n(M$_{\odot}$ pc$^2$ / yr)')


        if 'bh_pos_1' in plots:
            ax = next(na)
            ax.plot(time, d['BH_pos_1'][:, iaxis], label=r"${{{}}}_1$".format(axis))
            ax.set_ylabel('Pos 1')
        if 'bh_pos_2' in plots:
            ax = next(na)
            ax.plot(time, d['BH_pos_2'][:, iaxis], label=r"${{{}}}_2$".format(axis))
            ax.set_ylabel('Pos 2')

        if 'acc_torque' in plots:
            ax = next(na)
            ax.plot(time, d['acc_torque'][:, iaxis], label=r"$\tau_{{{}}}$".format(axis))
            ax.set_ylabel('Accretion Torque\n(M$_{\odot}$ pc$^2$ / yr$^{{-2}}$)')
            # ax.set_yscale('symlog', linthresh=1.e-13)
            ax.set_ylim([-5.e-12,5.e-12])
        if 'grav_torque' in plots:
            ax = next(na)
            ax.plot(time, d['grav_torque'][:, iaxis], label=r"$\tau_{{{}}}$".format(axis))
            ax.set_ylabel('Grav Torque\n(M$_{\odot}$ pc$^2$ / yr$^{{-2}}$)')
            # ax.set_yscale('symlog', linthresh=1.e-13)
            ax.set_ylim([-5.e-12,5.e-12])



        #         ax.plot(time,angmom[:,iaxis]/torque[:,iaxis],label=r"$'\tau_{{{}}}$".format(axis))
        #         ax.plot(time,d['net_force'][:,iaxis],label=r"$F_{{{}}}$".format(axis))
        if 'gcumtorque' in plots:
            ax = next(na)
            ax.plot(time, d['cum_torque'][:, iaxis], label=r"$\int\tau_{{{}}}dt$".format(axis))
            ax.set_ylabel('Integrated Grav Torque\n(M$_{\odot}$ pc$^2$ / yr)')


        if 'bh_vel_1' in plots:
            ax = next(na)
            ax.plot(time, d['BH_vel_1'][:, iaxis], label=r"$v_{{{},1}}$".format(axis))
            ax.set_ylabel('Vel 1')
        if 'bh_vel_2' in plots:
            ax = next(na)
            ax.plot(time, d['BH_vel_2'][:, iaxis], label=r"$v_{{{},2}}$".format(axis))
            ax.set_ylabel('Vel 2')

        if 'grav_force_1' in plots:
            ax = next(na)
            ax.plot(time, d['BH_force_1'][:, iaxis], label=r"$F_{{{},1}}$".format(axis))
            ax.set_ylabel('Grav Force 1')
        if 'grav_force_2' in plots:
            ax = next(na)
            ax.plot(time, d['BH_force_2'][:, iaxis], label=r"$F_{{{},2}}$".format(axis))
            ax.set_ylabel('Grav Force 2')

        if 'acc_force_1' in plots:
            ax = next(na)
            ax.plot(time, d['acc_force_1'][:, iaxis], label=r"$F_{{{},1}}$".format(axis))
            ax.set_ylabel('Acc Force 1')
        if 'acc_force_2' in plots:
            ax = next(na)
            ax.plot(time, d['acc_force_2'][:, iaxis], label=r"$F_{{{},2}}$".format(axis))
            ax.set_ylabel('Acc Force 2')

        if 'bh_circradforce_1' in plots:
            ax = next(na)
            ax.plot(time, (d['BH_force_1'][:, iaxis] - d['radial_force_1'][:, iaxis]),
                    label=r"$F_{{c,{},1}}$".format(axis))
            ax.plot(time, (d['radial_force_1'][:, iaxis]), label=r"$F_{{r,{},1}}$".format(axis))
            ax.set_ylabel('Radial/non-radial force 1')
        if 'bh_circradforce_2' in plots:
            ax = next(na)
            ax.plot(time, (d['BH_force_2'][:, iaxis] - d['radial_force_2'][:, iaxis]),
                    label=r"$F_{{c,{},2}}$".format(axis))
            ax.plot(time, (d['radial_force_2'][:, iaxis]), label=r"$F_{{r,{},2}}$".format(axis, axis))
            ax.set_ylabel('Radial/non-radial force 2')

        if 'bh_accmom' in plots:
            ax = next(na)
            ax.plot(time, d['BH_acc_mom_1'][:, iaxis], label=r"$P_{{a,{},1}}$".format(axis))
            ax.plot(time, d['BH_acc_mom_2'][:, iaxis], label=r"$P_{{a,{},2}}$".format(axis))
            ax.set_ylabel('Accreted momentum')

        if 'acc_cumtorque' in plots:
            ax = next(na)
            ax.plot(time, d['acc_torque'][:, iaxis], label=r"$\tau_{{{}}}$".format(axis))
            ax.set_ylabel('Accreted torque')

        if 'acc_cumtorque' in plots:
            ax = next(na)
            ax.plot(time, d['acc_cum_torque'][:, iaxis], label=r"$\int\tau_{{{}}}dt$".format(axis))
            ax.set_ylabel('Integrated Accretion Torque')

        if 'acc_J' in plots:
            ax = next(na)
            ax.plot(time, d['BH_acc_J_1'][:, iaxis]+d['BH_acc_J_2'][:, iaxis], label=r"$\int\tau_{{{}}}dt$".format(axis))
            ax.set_ylabel('Acc J (Msol pc$^{2}$ yr$^{-1}$)')

        if 'grav_J' in plots:
            ax = next(na)
            ax.plot(time, d['BH_grav_J_1'][:, iaxis]+d['BH_grav_J_2'][:, iaxis], label=r"$\int\tau_{{{}}}dt$".format(axis))
            ax.set_ylabel(r'Grav J (Msol pc$^{2}$ yr$^{-1}$)')

    #     sp[6,0].plot(time,(1.-d['BH_50_1']),label='1')
    #     sp[6,0].plot(time,(1.-d['BH_50_2']),label='2')

    #         sp[5,0].plot(time,d['BH_force_1'][:,iaxis]-d['BH_BH_accel_1'][:,iaxis],label=r"$\Delta F_{{{},1}}$".format(axis))
    #         sp[5,1].plot(time,d['BH_force_2'][:,iaxis]-d['BH_BH_accel_2'][:,iaxis],label=r"$\Delta F_{{{},2}}$".format(axis))
    #     sp[1,1].set_ylabel('Torque time-scale\n(yr)')
    #     sp[1,1].set_ylabel('Force sum')

    #     sp[6,0].set_ylabel('BH force 50 centile')
    #     sp[6,0].set_yscale('log')

    for ix in range(ncols):
        sp[-1, ix].set_xlabel(r'Time / yr')
        #         for iy in [4,5]:
        #             sp[iy,ix].set_yscale('symlog',linthresh=1.e-1)
        #             sp[iy,ix].set_ylim(-0.06,0.06)
        for iy in range(nrows):
            sp[iy, ix].legend()
            sp[iy, ix].set_xlim([0, tmax])

    fig.savefig(os.path.join(this_dir,f"../../figures/torque_{run_id}_{output_dir}.pdf"))
    plt.close('all')

