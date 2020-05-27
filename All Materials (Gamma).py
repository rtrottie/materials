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
from material_helpers import *
import numpy as np
import Make_Dos
try:
    client.close()
except:
    pass

(db,fs,client) = load_db()

#%% Import All Materials
# def get_matches(temp_match_criteria):
#     temp_matches = list(db.database.find(temp_match_criteria).sort('energy'))
#     key = 'defect_type'
#     if len(temp_matches) > 1:
#         view_multiple(temp_matches)
#         print('Multiple Matches {}'.format(temp_match_criteria))
#         [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
#         [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
#         print()
# #         view_multiple(temp_matches)
#         return temp_matches[0]
#     elif len(temp_matches) == 0:
#         print('No Matches {}'.format(temp_match_criteria))
#         return None
#     else:
#         return temp_matches[0]

match_criteria = {
    'pathway_count' : {'$exists' : True},
    'defect_type' : {'$exists' : False},
    'ts_type' : {'$exists' : False},
    'misc_labels' : 'kpts_div_2'

}

materials = [x['material'] for x in db.database.find(match_criteria).sort('material')]

print(len(materials))
# print(materials)


def get_matches(temp_match_criteria):
    temp_matches = list(db.database.find(temp_match_criteria).sort('energy'))
    key = 'defect_type'
    if len(temp_matches) > 1:
        #         view_multiple(temp_matches)
        print('Multiple Matches {}'.format(temp_match_criteria))
        [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        print()
        return temp_matches[0]
    elif len(temp_matches) == 0:
        print('No Matches {}'.format(temp_match_criteria))
        return None
    else:
        return temp_matches[0]


mat_dict = {}
ts_dict = {}

for material in materials:
    mat_dict[material[1]] = {}
    ts_dict[material[1]] = {}
    l = 7
    match_criteria = {
        'defect_location': {'$exists': False},
        'material': material,
        'labels': {'$nin': ['surface', 'ts', 'interpolation']},
        'pathway_count': {'$exists': True},
        'poscar.structure.lattice.a': {'$gt': l},
        'ts_type': {'$exists': False},
        'misc_labels' : 'kpts_div_2'
    }
    lb_criteria = {
        'material': material,
        'ts_type': 'pc',
        'index': -1,
        'misc_labels' : 'kpts_div_2'
    }
    ub_criteria = {
        'material': material,
        'interpolation_direction': {'$exists': True},
        'index': -1,
        'misc_labels' : 'kpts_div_2'
    }
    base = get_matches(match_criteria)

    match_criteria = {
        'defect_location': 'start',
        'material': material,
        'misc_labels' : 'kpts_div_2'
    }
    start = get_matches(match_criteria)

    mat_dict[material[1]]['base'] = base
    mat_dict[material[1]]['start'] = start
    mat_dict[material[1]]['finals'] = []
    ts_dict[material[1]]['lb'] = []
    ts_dict[material[1]]['ub'] = []
    for i in range(int(base['pathway_count'][0])):
        match_criteria['defect_location'] = 'final.{}'.format(i)
        lb_criteria['index'] = '{}'.format(i)
        ub_criteria['index'] = '{}'.format(i)

        final = get_matches(match_criteria)
        lb = get_matches(lb_criteria)
        ub = list(db.database.find(ub_criteria))

        mat_dict[material[1]]['finals'].append(final)
        ts_dict[material[1]]['lb'].append(lb)
        ts_dict[material[1]]['ub'].append(ub)



string = ''

fit = True

done = [
    'coga2o4',
    'cdal2o4',
    'laalo3',
    'ksbo3',
    'bacoo3',
    'bafeo3',
    'sralo3',
    'mgal2o4',
    'mnal2o4',
    'bazro3',
    'gemg2o4',
    'srfeo3',


]

bad_vacancies = [
    'bialo3',
]

converging = [
]

bad_energies = [
    'gefe2o4',
    'znfe2o4',
    'bialo3',
    'bigao3',
]

nup_wrong = [
]

checked = [
    'basno3',
    'ygao3',
    'vmg2o4',
    'cucr2o4',
    'inalo3',
    'alv2o4',
    'srmn3po3',
    'bimno3',
    'lagao3',
    'srsno3',
    'bageo3',
    'lamn3po3',
    'feal2o4',
    'cral2o4',
    'cual2o4',
    'yalo3',


]

get_both_sides = [
    'cucr2o4'
]

try2 = [
]

misaligned = [
    'znal2o4',
    'simg2o4',
    'sizn2o4',
    'feal2o4',
]

plot = [
]
aluminates = [
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

linewidth = 1.5
exact_range = 0.1

ranges = []
error_min = []
error_max = []
fig_all, ax_all = plt.subplots(figsize=(15, 15))
fig_pred, ax_pred = plt.subplots(figsize=(6, 5))
ax_pred.set_ylim([0, 4])
ax_pred.set_xlim([0, 4])

ax_all.plot([-1], [-1], color='blue', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.plot([-1], [-1], color='red', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.legend(['Perovskites', 'Spinels'], loc='upper left')
ax_all.set_ylabel('Diffusion Barrier (eV)')
ax_all.set_xlabel('Oxygen Vacancy Energy (eV)')

dont_plot = done + checked + try2 + misaligned + bad_energies + converging + nup_wrong

labels = [
    'laalo3',
    'lamn3po3',
    'sralo3',
    'srmn3po3',
]
plot = labels
plot = ['laalo3']

for material in mat_dict.keys():
    mat = mat_dict[material]
    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]
    plot_full(mat_dict, material, fig_all, ax_all, bars_list, labels=labels, label_all=label_all, string=string)
    if (plot_all or material in plot) and (material not in dont_plot):
        make_bar_chart(mat, bars_list, delta)


fig_all.show()
