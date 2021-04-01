import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# simulations = ["binary","unary"]
simulations = ["restest"]


with open("../local_config/skirt_data_dir") as f:
    skirt_data_dir = f.readline()

for sim in simulations:
    print(f"Summarising {sim}")
    # fname = f"{skirt_data_dir}/{sim}_lowres_allprobe_cellprops.dat"
    fname = f"{skirt_data_dir}/test_data/big_{sim}_allprobe_cellprops.dat"
    df = pd.read_csv(fname,sep=' ',header=None,skiprows=9,index_col=0,
                     names = ["x","y","z","V","tau","dust","e","H"])
    fname = f"{skirt_data_dir}/test_data/big_{sim}_tdust_T.dat"
    df_tdust = pd.read_csv(fname,sep=' ',header=None,skiprows=3,index_col=0,
                     names = ["Tdust"])
    df = pd.merge(df,df_tdust,left_index=True,right_index=True)

    # print(df)

    for key in df.keys():
        if np.sum(np.nonzero(df[key].values))==0:
            print(f"{key} is zero throughout run")
        else:
            print(f"{key}: Min:{np.nanmin(df[key].values)} Mean:{np.nanmean(df[key].values)} Median:{np.nanmedian(df[key].values)} Max:{np.nanmax(df[key].values)} ")
            print(f"{key}: Non-zero:{np.count_nonzero(df[key].values)/len(df)} Finite:{np.sum(np.isfinite(df[key].values))/len(df)}")

    print("Dust 2D histograms")

    df["r2d"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)
    for key in ["dust","Tdust","tau"]:
        # dust mean properties plots
        # H,xedges,yedges = np.histogram2d(df["z"],df["r2d"],weights=df[key],bins=1024)

        for coords,prefix in [[["y","x"],"yx"],[["z","r2d"],"rz"]]:

            H,yedges,xedges = np.histogram2d(df[coords[0]],df[coords[1]],weights=df[key],bins=128)
            fig,sp = plt.subplots(figsize=(6,6),dpi=200,constrained_layout=True)
            mp = sp.pcolormesh(xedges,yedges,H,norm=colors.LogNorm())
            fig.colorbar(mp)
            fig.savefig(f"../quick_out/{prefix}_{key}_{sim}.png")
            # fig.savefig(f"../quick_out/xy_{key}_{sim}.png")
            plt.close(fig)

