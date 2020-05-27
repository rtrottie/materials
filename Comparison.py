#%% Get Dictionaries
import matplotlib.pyplot as plt
exec(open(b'.\All Materials.py').read())
full_kpoints_thermo = mat_dict
full_kpoints_ts = ts_dict

exec(open(b'C:\Users\RyanTrottier\PycharmProjects\materials\All Materials (Gamma).py').read())
reduced_kpoints_thermo = mat_dict
reduced_kpoints_ts = ts_dict


#%% Get Info

redraw = False
fig_all, ax_all = plt.subplots(figsize=(10, 10))

ax_all.set_ylabel('Diffusion Barrier (eV)')
ax_all.set_xlabel('Oxygen Vacancy Energy (eV)')

common_materials = [ x for x in reduced_kpoints_thermo if x in full_kpoints_thermo ]

cmap = plt.cm.get_cmap('hsv', len(common_materials))
for i, material in enumerate(common_materials):
    mat = full_kpoints_thermo[material]
    thermo, lb_bars, ub_bars, delta = get_bars(mat, full_kpoints_ts[material])
    bars_list = [thermo, lb_bars, ub_bars]
    _,_, ts_estimate, ovac = plot_full(mat_dict, material, fig_all, ax_all, bars_list, colors=[cmap(i)]*3)


    mat = reduced_kpoints_thermo[material]
    thermo, lb_bars, ub_bars, delta = get_bars(mat, reduced_kpoints_ts[material])
    bars_list = [thermo, lb_bars, ub_bars]
    _, _, ts_estimate_red, ovac_red = plot_full(mat_dict, material, fig_all, ax_all, bars_list, colors=[cmap(i)]*3)
    if redraw:
        fig_all.show()
        input("Press Enter to continue...")
        fig_all, ax_all = plt.subplots(figsize=(5, 5))

        ax_all.set_ylabel('Diffusion Barrier (eV)')
        ax_all.set_xlabel('Oxygen Vacancy Energy (eV)')
    else:
        if ts_estimate < 99 and ts_estimate_red < 99:
            ax_all.plot([ovac, ovac_red], [ts_estimate, ts_estimate_red], linestyle=':', color=cmap(i) )

            ax_all.plot([ovac_red], [ts_estimate_red], marker='o', color=cmap(i))
            ax_all.plot([ovac], [ts_estimate], marker='x', color=cmap(i), markersize=8)

fig_all.show()