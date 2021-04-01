from generate_image_maps import t_images,t_orbit,runs,widths

import numpy as np
import h5py
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt


if __name__ == '__main__':
    cyberpunk_mode = True

    pic_inches = 3
    pic_dpi = 100

    if cyberpunk_mode:
        plt.style.use('dark_background')
        mpl.rcParams['font.family'] = 'monospace'
        cmap = "plasma"
    else:
        cmap = "cividis"

    f_in = h5py.File("../analysis_out/image_maps_evolution.hdf5", "r")

    # all_vmin = None
    # all_vmax = None
    # for (run_dir, run_names),t_image in zip(runs,t_images):
    #     for run_name in run_names:
    #         for iy in range(len(t_image)):
    #             key = f"{run_dir}_{run_name}_{iy}"
    #             m = np.array(f_in[key])
    #             m = m[np.isfinite(m)]
    #             if all_vmin is None:
    #                 all_vmin = np.nanmin(m)
    #                 all_vmax = np.nanmax(m)
    #             else:
    #                 all_vmin = min(np.nanmin(m),all_vmin)
    #                 all_vmax = max(np.nanmax(m),all_vmax)
    #
    # print(all_vmin,all_vmax)

    # all_vmin = 1.e-25
    # all_vmax = 1.e-15
    all_vmin = -25
    all_vmax = -15


    def plot_key(ax, key, width) :
        extent = [-width / 2, width / 2, -width / 2, width / 2]
        m = np.log10(np.array(f_in[key]))
        # print(key)
        # print(np.nanmax(m), np.nanmin(m[np.isfinite(m)]))
        ax.imshow(m, vmin=all_vmin, vmax=all_vmax, cmap=cmap, extent=extent)
        # ax.imshow(m, norm=colors.LogNorm(vmin=all_vmin, vmax=all_vmax), cmap=cmap, extent=extent)

    # x = [1,2,3]
    # y = [3,4,5]
    #
    # x2 = [5,6,7]
    # y2 = [6,5,4]

    for i_fig,((run_dir, run_names),t_image) in enumerate(zip(runs,t_images)):
        nx = len(run_names)*len(widths)*2 # 2 for face/edge
        ny = len(t_image)
        # fig,sp = plt.subplots(ny,nx,figsize=(pic_inches*nx,pic_inches*ny),dpi=pic_dpi,constrained_layout=True,sharex='col',sharey='col')
        fig,sp = plt.subplots(ny,nx,figsize=(pic_inches*nx,pic_inches*ny),dpi=pic_dpi,constrained_layout=True,sharex='col',sharey='col')
        for irun,run_name in enumerate(run_names):
            for iy in range(len(t_image)):
                print(f"Plotting {run_dir} {run_name} {iy}")
                for iwidth,width in enumerate(widths):
                    key = f"{run_dir}_{run_name}_{iy}_{iwidth}_face"
                    ix = irun*4 + iwidth
                    print(ix,irun,iwidth,"face")
                    ax = sp[iy,ix]
                    plot_key(ax,key,width)
                    # ax.scatter(x,y)

                    key = f"{run_dir}_{run_name}_{iy}_{iwidth}_edge"
                    ix = irun*4 + 2 + iwidth
                    print(ix,irun,iwidth,"edge")
                    ax = sp[iy,ix]
                    plot_key(ax,key,width)
                    # ax.scatter(x2, y2)
                # ax = sp[iy,-1]
                # # ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                # # ax2 = ax.twinx()
                ax.set_ylabel(f't={t_image[iy]}T')
        print("Saving")
        fig.savefig(f"../figures_out/evolution_{i_fig}.png")
        plt.close(fig)
    f_in.close()
