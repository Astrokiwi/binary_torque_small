from generate_image_maps import t_images,t_orbit,runs,widths,major_axis
from gtools import gizmo_tools

import numpy as np
import h5py
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import colors
import matplotlib.pyplot as plt
import copy


if __name__ == '__main__':
    pic_inches = 1.5
    pic_dpi = 200

    cyberpunk_mode = False

    box_color='green'
    if cyberpunk_mode:
        plt.style.use('dark_background')
        mpl.rcParams['font.family'] = 'monospace'
        cmap = "plasma"
    else:
        cmap = "plasma"
        # cmap = "cividis"

    current_cmap = copy.copy(mpl.cm.get_cmap(cmap))
    current_cmap.set_bad(color='black')

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
        ax.imshow(m, vmin=all_vmin, vmax=all_vmax, cmap=current_cmap, extent=extent)
        # ax.imshow(m, norm=colors.LogNorm(vmin=all_vmin, vmax=all_vmax), cmap=cmap, extent=extent)

    # x = [1,2,3]
    # y = [3,4,5]
    #
    # x2 = [5,6,7]
    # y2 = [6,5,4]

    for i_fig,((run_dir, run_names,tidy_names),t_image) in enumerate(zip(runs,t_images)):
        for irun,(run_name,tidy_name) in enumerate(zip(run_names,tidy_names)):
            nx = len(widths) * 2
            ny = len(t_image)
            # fig,sp = plt.subplots(ny,nx,figsize=(pic_inches*nx,pic_inches*ny),dpi=pic_dpi,constrained_layout=True,sharex='col',sharey='col')
            fig, sp = plt.subplots(ny, nx, figsize=(pic_inches * nx, pic_inches * ny), dpi=pic_dpi,
                                   constrained_layout=True, sharex='col', sharey='col')
            for iy in range(len(t_image)):
                print(f"Plotting {run_dir} {run_name} {iy}")
                for iwidth,width in enumerate(widths):
                    scaled_width = width/major_axis
                    key = f"{run_dir}_{run_name}_{iy}_{iwidth}_face"
                    ix = iwidth
                    print(ix,irun,iwidth,"face")
                    ax = sp[iy,ix]
                    plot_key(ax,key,scaled_width)
                    if ix!=0:
                        ax.set_yticklabels([])
                    # ax.scatter(x,y)

                    key = f"{run_dir}_{run_name}_{iy}_{iwidth}_edge"
                    ix = 2 + iwidth
                    print(ix,irun,iwidth,"edge")
                    ax = sp[iy,ix]
                    plot_key(ax,key,scaled_width)
                    ax.set_yticklabels([])
                    # ax.scatter(x2, y2)


                # ax = sp[iy,-1]
                # # ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                # # ax2 = ax.twinx()
                ax.set_ylabel(f't={t_image[iy]}T')
            box_corner = [-widths[0]*.45/major_axis,-widths[0]*.45/major_axis]
            box_dims = [widths[0]*.9/major_axis,widths[0]*.9/major_axis]
            for iy in range(len(t_image)):
                gizmo_tools.box_connected_two_axes(sp[iy, 0], sp[iy, 1], box_corner, box_dims,color=box_color,corners=[[1,0],[1,1]])
                gizmo_tools.box_connected_two_axes(sp[iy, 2], sp[iy, 3], box_corner, box_dims,color=box_color,corners=[[1,0],[1,1]])

            fig.suptitle(tidy_name)
            print("Saving")

            fig.savefig(f"../figures_out/evolution_{run_dir}_{run_name}.png")
            plt.close(fig)
    f_in.close()
