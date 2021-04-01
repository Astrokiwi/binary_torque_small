from gtools import minidisc_masses

if __name__ == '__main__':
    runs = [["rad_prod",["rad_small_ecc_earlier","rad_small_circ_earlier"]],
            ["rad_prod",["rad_small_ecc_full","rad_small_circ_full"]],
            ["norad_prod",["small_ecc","small_circ"]]]
    major_axis = 0.035

    for directory,sims in runs:
        for sim in sims:
            minidisc_masses.analyse_run(directory,sim,major_axis,"../analysis_out")