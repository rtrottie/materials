#%% Import

import matplotlib
import matplotlib.pyplot as plt
import pymongo
import gridfs
import bz2
import Vis
import time
from bson import ObjectId
from Classes_Pymatgen import *
from AddDB import *
from View_Structures import *
from Database_Tools import get_file
from Helpers import get_FERE_chemical_potential
from scipy.stats import gmean
import numpy as np
import Make_Dos
import copy
from material_helpers import *
try:
    client.close()
except:
    pass

(db,fs,client) = load_db()

#%% Import All Materials

verbose = False

match_criteria = {
    'pathway_count' : {'$exists' : True},
    'defect_type' : {'$exists' : False},
    'ts_type' : {'$exists' : False}
}

materials = [x['material'] for x in db.database.find(match_criteria).sort('material')]
# materials = [
# ['perovskite', 'bazro3'],
#              ['spinel', 'mnal2o4']
# ]
if verbose:
    print(len(materials))
# print(materials)


def get_matches(temp_match_criteria):
    temp_matches = list(db.database.find(temp_match_criteria).sort('energy'))
    key = 'defect_type'
    if len(temp_matches) > 1:
        #         view_multiple(temp_matches)
        # print('Multiple Matches {}'.format(temp_match_criteria))
        # [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        # [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        # print()
        return temp_matches[0]
    elif len(temp_matches) == 0:
        # print('No Matches {}'.format(temp_match_criteria))
        return None
    else:
        return temp_matches[0]


mat_dict = {}
ts_dict = {}

for material in materials:
    l = 7
    match_criteria = {
        'defect_location': {'$exists': False},
        'material': material,
        'labels': {'$nin': ['surface', 'ts', 'interpolation', 'charged_defect']},
        'pathway_count': {'$exists': True},
        'poscar.structure.lattice.a': {'$gt': l},
        'ts_type': {'$exists': False},
        'energy' : {'$exists' : True},
        'misc_labels' : {'$nin' : ['kpts_div_2']},
    }
    base = get_matches(match_criteria)
    lb_criteria = {
        'labels': {'$nin': ['surface', 'interpolation'], '$all' : ['charged_defect']},
        'material': material,
        'ts_type': 'pc',
        'index': -1,
        'misc_labels' : {'$nin' : ['kpts_div_2']},
        'energy' : {'$exists' : True},
        'kpoints.kpoints' : base['kpoints']['kpoints'],
    }
    ub_criteria = {
        'labels': {'$nin': ['surface'], '$all' : ['charged_defect']},
        'material': material,
        'interpolation_direction': {'$exists': True},
        'index': -1,
        'misc_labels' : {'$nin' : ['kpts_div_2']},
        'kpoints.kpoints' : base['kpoints']['kpoints'],
    }

    match_criteria = {
        'labels': {'$nin': ['surface', 'ts', 'interpolation'], '$all' : ['charged_defect']},
        'defect_location': 'start',
        'material': material,
        'misc_labels' : {'$nin' : ['kpts_div_2']},
        'energy' : {'$exists' : True},
        'kpoints.kpoints' : base['kpoints']['kpoints'],
    }
    start = get_matches(match_criteria)

    if start:

        lb0_criteria = {
            'material': material,
            'ts_type': 'pc',
            'index': -1,
            'misc_labels': {'$nin': ['kpts_div_2']},
            'energy': {'$exists': True},
            'kpoints.kpoints': base['kpoints']['kpoints'],
            'labels': {'$nin': ['surface', 'interpolation', 'charged_defect']},
        }
        ub0_criteria = {
            'material': material,
            'interpolation_direction': {'$exists': True},
            'index': -1,
            'misc_labels': {'$nin': ['kpts_div_2']},
            'kpoints.kpoints': base['kpoints']['kpoints'],
            'labels': {'$nin': ['surface', 'charged_defect']},
        }

        match_criteria0 = {
            'defect_location': 'start',
            'material': material,
            'misc_labels': {'$nin': ['kpts_div_2']},
            'energy': {'$exists': True},
            'kpoints.kpoints': base['kpoints']['kpoints'],
            'labels': {'$nin': ['surface', 'ts', 'interpolation', 'charged_defect']},
        }
        start0 = get_matches(match_criteria0)

        mat_dict[material[1]] = {}
        ts_dict[material[1]] = {}

        mat_dict[material[1]]['base'] = base
        mat_dict[material[1]]['start'] = start
        mat_dict[material[1]]['start0'] = start0
        mat_dict[material[1]]['finals'] = []
        ts_dict[material[1]]['lb'] = []
        ts_dict[material[1]]['ub'] = []
        mat_dict[material[1]]['finals0'] = []
        ts_dict[material[1]]['lb0'] = []
        ts_dict[material[1]]['ub0'] = []
        for i in range(int(base['pathway_count'][0])):
            match_criteria['defect_location'] = 'final.{}'.format(i)
            match_criteria0['defect_location'] = 'final.{}'.format(i)
            lb_criteria['index'] = '{}'.format(i)
            ub_criteria['index'] = '{}'.format(i)
            lb0_criteria['index'] = '{}'.format(i)
            ub0_criteria['index'] = '{}'.format(i)

            final = get_matches(match_criteria)
            final0 = get_matches(match_criteria0)
            lb = get_matches(lb_criteria)
            ub = list(db.database.find(ub_criteria))
            lb0 = get_matches(lb0_criteria)
            ub0 = list(db.database.find(ub0_criteria))

            mat_dict[material[1]]['finals'].append(final)
            mat_dict[material[1]]['finals0'].append(final0)
            ts_dict[material[1]]['lb'].append(lb)
            ts_dict[material[1]]['ub'].append(ub)
            ts_dict[material[1]]['lb0'].append(lb0)
            ts_dict[material[1]]['ub0'].append(ub0)




#%% Plotting Functions

string = ''

fit = True

done = [

]

done_lower_ub = [
]

converging = [
]

bad_energies = [
]

nup_wrong = [
]

lookatit = [
]


checked = [
]

try2 = [
]

misaligned = [
    'bisbo3',
]

no_nupdown = [
]
gone = [
]
plot = [
]
aluminates = [
]
no_vacanciestxt = [
]



#%% Plot All

plot_all = True
label_all = False

try:
    for fig in figs:
        plt.close(fig)
    print('Cleaned Figs')
except:
    figs = []

linewidth = .75
exact_range = 0.1

ranges = []
error_min = []
error_max = []
fig_all, ax_all = plt.subplots(figsize=(3, 3))
fig_pred, ax_pred = plt.subplots(figsize=(6, 5))
ax_pred.set_ylim([0, 4])
ax_pred.set_xlim([0, 4])

ax_all.plot([-1], [-1], color='blue', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.plot([-1], [-1], color='red', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.legend(['Perovskites', 'Spinels'], loc='upper left')
ax_all.set_ylabel('Diffusion Barrier (eV)')
ax_all.set_xlabel('Oxygen Vacancy Energy (eV)')

ax_all.set_ylim([0,4])
ax_all.set_xlim([-2, 10])

dont_plot = done + done_lower_ub + checked + try2 + misaligned + bad_energies + converging + nup_wrong + no_nupdown + no_vacanciestxt + gone
plot = ['bialo3']

labels = [
    'laalo3',
    'lamn3po3',
    'sralo3',
    'srmn3po3',
]
labels = []
csv_dict = {}

csv_dict = Csv_dict()
label_list = []
lb = []
lb0 = []
avgs = []
avgs0 = []
err = []
err0 = []
minmax = 0.05

unit_criteria = {
            'defect_type': {'$exists': False},
            'ts_type': {'$exists': False},
            'material': {'$all': [], '$in': ['spinel', 'perovskite'], '$nin': ['from_zach']},
            'labels' : {
                '$all' : ['unit_cell'],
                '$nin' : ['surface']}
        }

if __name__ == '__main__':

    for material in mat_dict.keys():
    # for material in ['cual2o4']:
        mat = mat_dict[material]
        thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
        thermo0, lb0_bars, ub0_bars, delta = get_bars(mat, ts_dict[material], suffix='0')
        if material == 'bazro3':
            pass
        if verbose:
            print(lb_bars, lb0_bars)
        bars_list = [thermo, lb_bars, ub_bars]


        ranges = np.array(ub_bars) - np.array(lb_bars)

        ranges = np.array([ x if x > 0 else 0 for x in ranges])
        avg = np.array(lb_bars) + ranges/2 #type:np.array
        if max(avg) > 0:
            label_list.append(material)
            ranges0 = np.array(ub0_bars) - np.array(lb0_bars)
            ranges0 = np.array([ x if x > 0 else 0 for x in ranges0])
            avg0 = np.array(lb0_bars) + ranges0/2
            i = list(avg).index(min([ x for x in avg if x > 0 ]))
            lb.append(lb_bars[i])
            err.append(max(ranges[i], minmax))
            i0 = list(avg0).index(min([ x for x in avg0 if x > 0 ]))
            lb0.append(lb0_bars[i0])
            err0.append(max(ranges0[i0], minmax))
            avgs.append(avg[i0])
            avgs0.append(avg0[i0])
        elif min(avg) < 0:
            pass
            # del mat_dict[material]
            # del ts_dict[material]
        else:
            label_list.append(material)
            lb.append(0)
            lb0.append(0)
            err.append(0)
            err0.append(0)

        if (err[-1] == 0 or lb[-1] == 0) and (material not in misaligned):
            print(label_list[-1])
            ranges0 = np.array(ub0_bars) - np.array(lb0_bars)
            ranges0 = np.array([ x if x > 0 else 0 for x in ranges0])
            avg0 = np.array(lb0_bars) + ranges0 / 2
            print(list(avg0).index(min([ x for x in avg0 if x > 0 ])))
            print(err[-1])
            print(lb[-1])

    figsize = (10,5)
    font = {
        'family' : 'Arial',
        'size' : 12,
    }
    label_font = copy.deepcopy(font)
    label_font['size'] = 14

    tick_font = copy.deepcopy(font)
    tick_font['size'] = 9
    labels = {
        'bacoo3' : 'BaCoO$_3$',
        'bafeo3' : 'BaFeO$_3$'
    }
    def make_plot(bars_list, bottoms, legend=[], colors=['k', 'r', 'b'], ylim=None, ylabel='',label_list=label_list):
        width = .15 / len(bars_list)
        x = np.arange(len(bars_list[0]))
        fig, ax = plt.subplots(figsize=figsize)
        fig.patch.set_facecolor('white')

        for i, bars in enumerate(bars_list):
            rects = plt.errorbar(x + i * width, bars, color=colors[i], label=legend[i], yerr=np.array(bottoms[i])/2,
                                 fmt='o', markersize=0, linewidth=10)
        lables = ax.set_xticklabels([ labels[x] if x in labels else x for x in label_list], fontdict=tick_font)
        ax.set_ylabel(ylabel, fontdict=label_font)

        ticks = ax.set_xticks(x + width * (len(bars_list)-1) / 2)
        # ticks = ax.set_xticks(x)
        if ylim:
            ax.set_ylim(ylim)
        plt.show()

    make_plot([lb, lb0], [err, err0], legend=['+2', '0'])

    #%%

    maxrange = 0.5
    plt.subplots(figsize=(3,3))
    x = np.array(lb0) + np.array(err0)/2
    y = np.array(lb) + np.array(err)/2
    xerr = np.array(err0)/2
    yerr = np.array(err)/2

    x = x[np.where(xerr < maxrange)]
    y = y[np.where(xerr < maxrange)]
    yerr = yerr[np.where(xerr < maxrange)]
    ll = np.array(label_list)[np.where(xerr < maxrange)]

    xerr = xerr[np.where(xerr < maxrange)]




    color = []
    for material in ll:
        print(material)
        f=get_file(fs, mat_dict[material]['base']['vasprun'])
        vasprun = Vasprun(f)
        color.append(vasprun.get_band_structure().get_band_gap()['energy'])
        os.remove(f)
        # color.append(mat_dict[material]['base']['descriptors']['bandgap'])

    # #%%
    #
    maxbarrier = 5.5
    plt.plot([0, maxbarrier], [0, maxbarrier], '--', color='gray')
    plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='none', color='k')
    scat = plt.scatter(x, y, c=color, cmap='gnuplot2', marker='o')
    # scat = plt.scatter(x, y, marker='o)
    plt.colorbar(scat, label='Band Gap (eV)')
    # plt.colorbar(scat, label='V$_O$ Energy (eV)')
    plt.xlabel('V$_O$ Diffusion Barrier (eV)')
    plt.ylabel('V$_O^{2+}$ Diffusion Barrier (eV)')
    plt.xlim([0,maxbarrier])
    plt.ylim([0,maxbarrier])
    # for i in range(len(ll)):
    #     plt.text(x[i], y[i], ll[i])


    plt.show()
    #%%

    # for i in range(len(err)):
    #     print(i)
    #     if err[i] == 0 or lb[i] ==0:
    #         print(label_list[i])
    #         # print(i)
    #         print(err[i])
    #         print(lb[i])