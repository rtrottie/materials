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
try:
    client.close()
except:
    pass

(db,fs,client) = load_db()

#%% Import All Materials
def get_matches(temp_match_criteria):
    temp_matches = list(db.database.find(temp_match_criteria).sort('energy'))
    key = 'defect_type'
    if len(temp_matches) > 1:
        view_multiple(temp_matches)
        print('Multiple Matches {}'.format(temp_match_criteria))
        [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        print()
#         view_multiple(temp_matches)
        return temp_matches[0]
    elif len(temp_matches) == 0:
        print('No Matches {}'.format(temp_match_criteria))
        return None
    else:
        return temp_matches[0]

match_criteria = {
    'pathway_count' : {'$exists' : True},
    'defect_type' : {'$exists' : False},
    'ts_type' : {'$exists' : False}
}

materials = [x['material'] for x in db.database.find(match_criteria).sort('material')]

print(len(materials))
print(materials)


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
        'ts_type': {'$exists': False}
    }
    lb_criteria = {
        'material': material,
        'ts_type': 'pc',
        'index': -1
    }
    ub_criteria = {
        'material': material,
        'interpolation_direction': {'$exists': True},
        'index': -1
    }
    base = get_matches(match_criteria)

    match_criteria = {
        'defect_location': 'start',
        'material': material
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


#%% Plot
def get_band_center(vasprun, orbital, element='O'):
    (title, columns, scaling_factors) = Make_Dos.make_dos(vasprun, ['{}:{}'.format(orbital, element)])
    maximum = np.argmax(np.array(columns[0]) > 0)
    energies = columns[0]
    dos = np.array(columns[3]) - np.array(columns[4])
    return np.average(energies[:maximum], weights=dos[:maximum])


#     if not elements:
#         elements = vasprun.final_structure.species
#     maximum = np.argmax(vasprun.tdos.energies>0)
#     energies = vasprun.tdos.energies[:maximum]
#     centers = []
#     weights = []
#     dos = vasprun.complete_dos
#     total_dos = np.array([])
#     for i, a in enumerate( [ x for x in vasprun.final_structure if (x.specie in elements)]):
#         for orbital in orbitals:
#             site_dos = Make_Dos.get_dos(vasprun.complete_dos, i, 'p')
#             if len(total_dos) > 0:
#                 total_dos += site_dos.get_densities()
#             else:
#                 total_dos = site_dos.get_densities()

#     return np.average(site_dos.energies[:maximum], weights=total_dos[:maximum])

def add_to_string(string, base, ovac, bounds, weight=2, max_val=2):
    labels = ('material', 'barrier', 'ovac', 'density', 'fermi', 'o_p_band')
    labels = ('material', 'barrier', 'density', 'fermi', 'o_p_band')
    material = base['material'][1]
    if 'descriptors' in base and all([x in base['descriptors'] for x in labels[3:]]):  # if all labels are made
        descriptors = base['descriptors']
    else:
        s = Structure.from_dict(base['poscar']['structure'])
        v_file = get_file(fs, base['vasprun'])
        v = Vasprun(v_file)
        os.remove(v_file)
        # set descriptors
        descriptors = {
            'density': s.volume / len(s),
            'bandgap': v.tdos.get_gap(),
            'fermi': v.efermi,
            'o_p_band': get_band_center(v, 'p', 'O'),
            'p_band': get_band_center(v, 'p', 'O')
        }
        base['descriptors'] = descriptors
        db.database.update_one({'_id': base['_id']}, {'$set': {'descriptors': descriptors}})
    properties = [descriptors[label] for label in labels[2:]]

    predicted = (
            0.73118E+00 +
            0.96547E-04 * ovac ** 6 +
            0.75320E-02 * np.exp(ovac) * ovac / descriptors['o_p_band']

        #         -0.28972E+01 +
        #         0.45890E+02 * np.exp(descriptors['o_p_band']) * descriptors['density'] *  descriptors['density']

    )

    format_string = '{}  ' + '{:.4f}  ' * len(labels[1:]) + '\n'

    if string == '':
        label_string = '{}   ' * len(labels) + '\n'
        string = label_string.format(*labels)

    difference = abs(bounds[0] - bounds[1])
    if difference > max_val:
        return string, predicted
    else:
        weight = max(2, weight * (max_val - difference) / max_val)
    for bound in np.linspace(bounds[0], bounds[1], weight):
        string += format_string.format(material, bound, *properties)
    #         string += format_string.format(material, bound, ovac, *properties)
    return string, predicted


string = ''

fit = True

done = [
    'bacoo3',
    'bafeo3',
    'bamn3po3',
    'basno3',
    'bazro3',
    'bimno3',
    'kbio3',
    'ksbo3',
    'lagao3',
    'lamn3po3',
    'mgcr2o4',
    'sralo3',
    'srcoo3',
    'srfeo3',
    'srsno3',
    'ygao3',

    'cdal2o4',
    'coal2o4',
    'coga2o4',
    'cral2o4',
    'cual2o4',
    'feal2o4',
    'mgal2o4',
    'mnal2o4',
    'nial2o4',
    'vmg2o4',

]

converging = [
    'cafeo3',
    'ceo2',

    'cucr2o4',
    'crcr2o4',
]

bad_energies = [
    'tife2o4',
    'nife2o4',
    'mnv2o4',
    'tico2o4',
    'fecr2o4',
    'mnmn2o4',
    'fev2o4',
    'vfe2o4',
    'nini2o4',

    'snbeo3',
    'yfeo3',
    'bicoo3',
    'zrbeo3',
    'tibeo3',
]

nup_wrong = [
    'fega2o4'
]

checked = [
    'bageo3',
    'baseo3',
    'bisco3',
    'bivo3',
    'cageo3',
    'inalo3',
    'laalo3',
    #     'lagao3',
    'nasbo3',
    'snfeo3',
    'srmn3po3',
    'yalo3',

    'alal2o4',
    'alv2o4',
    'femn2o4',
    #     'gefe2o4',
]

try2 = [
]

misaligned = [
    'canio3',
    'cocr2o4',
    'geni2o4',
    'bisbo3',
    'bigao3',
    'caseo3',
    'limn2o4',
    'cuga2o4',
    'bialo3',
    'simg2o4',
    'gemg2o4',
    'sizn2o4',
    'sicd2o4',
    'navo3',
    'znal2o4',
]

plot = [
    'snfeo3',
    'srfeo3',
    'srcoo3',
    'bafeo3',
    'bacoo3',
    'lafeo3',
    'lacoo3',
]
aluminates = [
    'mnal2o4',
    'coal2o4',
    'nial2o4'
]

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
fig_all, ax_all = plt.subplots(figsize=(6, 5))
fig_pred, ax_pred = plt.subplots(figsize=(6, 5))
ax_pred.set_ylim([0, 4])
ax_pred.set_xlim([0, 4])

ax_all.plot([-1], [-1], color='blue', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.plot([-1], [-1], color='red', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.legend(['Perovskites', 'Spinels'], loc='upper left')
ax_all.set_ylabel('Diffusion Barrier (eV)')
ax_all.set_xlabel('Oxygen Vacancy Energy (eV)')

dont_plot = done + checked + try2 + misaligned + bad_energies + converging + nup_wrong

for material in mat_dict.keys():
    #   if material in aluminates:
    #   if material not in dont_plot:
    if material in done:
        mat = mat_dict[material]
        lb_list = ts_dict[material]['lb']
        ub_list = ts_dict[material]['ub']
        width = 0.75

        bulk_e = mat['base']['energy']

        bars = [mat['start']['energy'] if type(mat['start']) == dict else 0] + [x['energy'] if type(x) == dict else 0
                                                                                for x in mat['finals']]

        lb_bars = [0]
        delta = np.array([0])
        for i, energy in enumerate(bars):
            if i == 0:
                start_energy = energy
            else:
                lb_bars.append(
                    lb_list[i - 1]['energy'] - min(energy, start_energy) if energy < 0 and start_energy < 0 and type(
                        lb_list[i - 1]) == dict else 0)
                delta = np.append(delta, np.array(abs(start_energy - energy)))

        ub_bars = [0]
        for i, energy in enumerate(bars):
            if i == 0:
                start_energy = energy
            else:
                ub_bars.append(max([max(x['energies']) for x in ub_list[i - 1]]) - min(energy,
                                                                                       start_energy) if energy < 0 and start_energy < 0 and len(
                    ub_list[i - 1]) > 0 and type(ub_list[i - 1][0]) == dict else 0)

        bars = [(x + get_FERE_chemical_potential('O')) - bulk_e if x != 0 else 0 for x in bars]
        bars_list = [bars, lb_bars, ub_bars]
        colors = ['red', 'blue', 'green', 'yellow']
        xlabels = ['Start'] + ['Final.{}'.format(i) for i in range(len(bars) - 1)]
        xlabels = ['Start'] + [
            '{} through {}'.format(x['pathway_id']['initial-final'], x['pathway_id']['critical-triangle']) if type(
                x) == dict and 'pathway_id' in x else 'Final.{}'.format(i) for i, x in enumerate(mat['finals'])]
        legend = ['Ovac', 'LB Diffusion', 'UB Diffusion']
        width = .75 / len(bars_list)
        x = np.arange(len(bars))
        fig, ax = plt.subplots(figsize=(10, 3))
        for i, bars in enumerate(bars_list):
            if i >= 1:
                lower = np.array(bars) - delta
                lower = [bars[i] - delta[i] if bars[i] > 0 else 0 for i in range(len(bars))]
                delta_temp = [delta[i] if bars[i] > 0 else 0 for i in range(len(bars))]
                rects = ax.bar(x + i * width, lower, width, color=colors[i])
                rects = ax.bar(x + i * width, delta_temp, width, color='dark{}'.format(colors[i]), bottom=lower,
                               label='_nolegend_')
            else:
                rects = ax.bar(x + i * width, bars, width, color=colors[i])

        ts_estimates = [gmean([bars_list[1][i], bars_list[2][i]]) if bars_list[2][i] else 9999 for i in
                        range(len(bars_list[0]))]
        lowest_ts_estimate = min(ts_estimates)
        for i in range(len(bars_list[0])):
            if bars_list[0][i] and bars_list[1][i] and bars_list[2][i]:
                if 'spinel' in mat_dict[material]['base']['material']:
                    color = 'red'
                elif 'perovskite' in mat_dict[material]['base']['material']:
                    color = 'blue'
                zorder = 2
                if ts_estimates[i] > lowest_ts_estimate:
                    zorder = 1
                    if color == 'red':
                        color = '#F5B7B1'
                    elif color == 'blue':
                        color = '#AED6F1'
                    else:
                        color = 'light' + color
                else:
                    ranges.append(-bars_list[1][i] + bars_list[2][i])
                    if fit:
                        string, predicted = add_to_string(string, mat_dict[material]['base'], bars_list[0][i],
                                                          (bars_list[1][i], bars_list[2][i]),
                                                          weight=3, max_val=1.5)
                        error_min.append(min(abs(bars_list[1][i] - predicted), abs(predicted - bars_list[2][i])))
                        error_max.append(max(abs(bars_list[1][i] - predicted), abs(predicted - bars_list[2][i])))
                        ax_pred.plot([predicted, predicted], [bars_list[1][i], bars_list[2][i]], marker='_',
                                     color=color,
                                     linewidth=linewidth, markersize=8, markeredgewidth=linewidth, zorder=zorder)
                #                     ax_all.plot([bars_list[0][i]], [ predicted ], marker='o', color=color, markeredgecolor=color, markerfacecolor='none',
                #                              linewidth=linewidth, markersize=5, markeredgewidth=linewidth, zorder=zorder)

                ax_all.plot([bars_list[0][i], bars_list[0][i]], [bars_list[1][i], bars_list[2][i]], marker='_',
                            color=color,
                            linewidth=linewidth, markersize=8, markeredgewidth=linewidth, zorder=zorder)

        if 'other_work' in mat['base']:
            for i, migration_barrier in enumerate(mat['base']['other_work']['migration_barriers']):
                if 'energy' in migration_barrier:
                    energy = migration_barrier['energy']
                    label = '{} through {}'.format(migration_barrier['initial-final'],
                                                   migration_barrier['critical-triangle'])
                    if label in xlabels:
                        index = i + 1
                        ax.plot([index + width, index + width * len(bars_list)], [energy, energy], linewidth=2,
                                color='teal', label='_nolegend_')

        lables = ax.set_xticklabels(range(len(bars)))
        ax.set_ylabel('Energy (eV)')
        ax.set_title(material)
        ax.legend(legend, bbox_to_anchor=(0.5, -0.1), loc='upper center')
        ticks = ax.set_xticks(x + width * len(bars_list) / 2)
        ax.set_xticklabels(xlabels)
        time.sleep(0.001)
        figs.append(fig)

#     x = np.arange(max_i)
#     fig, ax = plt.subplots(figsize=(8,3))
#     for i, bars in enumerate(bars_list):
#         rects = ax.bar(x+i*width, bars, width, color=colors[i])
#     lables = ax.set_xticklabels(range(max_i))
#     ax.legend(legend,bbox_to_anchor=(0.5,-0.1),loc='upper center')


# /scratch/rtrottie/vasp/watersplitting/materials/10/pbeu/bulk/perovskite/baseo3/pc.5/0/neb.2