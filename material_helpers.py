import numpy as np
from Helpers import get_FERE_chemical_potential
from Classes_Pymatgen import *
import os
from scipy.stats import gmean
import matplotlib.pyplot as plt
from Database_Tools import get_file
import csv

font_dict = {'family' : 'Arial', 'size' : 16}
label_font = font_dict.copy()
tick_font = font_dict.copy()
tick_font['size'] = tick_font['size'] - 4
def get_lower_material(run):
    lower = run['material'][1]
    if lower == 'spinel' or lower == 'perovskite':
        lower = run['material'][0]
    return lower

class Csv_dict():

    def __init__(self):
        return

    materials = {}

    def get_property(self, run, prop):
        return run[prop]

    def get_energy(self, run, start, final):
        return run['energy']

    def get_material(self, run, start=None, final=None):
        lower = get_lower_material(run)
        lower = lower.replace('3p', '')
        lower = lower.replace('o3', 'O3')
        lower = lower.replace('o4', 'O4')
        indices = []
        elements = run['elements'].copy()
        elements.pop(elements.index('O'))
        for element in elements:
            indices.append(lower.find(element.lower()))
        next_indices = [indices[i]+len(elements[i])-1 for i in range(len(indices)) if len(elements[i]) > 1]
        # print(elements)
        # print(lower)
        # print(indices)
        # print(next_indices)
        for i in next_indices:
            if i in indices:
                raise Exception('Failed')
        for element in elements:
            lower = lower.replace(element.lower(), element)
        # print(lower)
        return lower

    def get_type(self, run, start, final):
        return run['material'][0]

    properties = [
        ('type',  get_type),
        # ('energy', get_energy)
    ]

    ts_properties = [

    ]

    def get_line(self, material, i):
        line = []
        line.append(material)
        line.append(i)
        line.append(self.materials[material][i]['ovac'])
        line.append(self.materials[material][i]['lb'])
        line.append(self.materials[material][i]['ub'])
        for property, _ in self.properties:
            line.append(self.materials[material][i][property])
        return line

    def get_header(self):
        line = ['material', 'i', 'ovac', 'lb', 'ub']
        for property, _ in self.properties:
            line.append(property)
        return line

    def get_list(self):
        lines = []
        lines.append(self.get_header())
        for material in self.materials.keys():
            for i in self.materials[material].keys():
                lines.append(self.get_line(material, i))
        return lines

    def get_column(self, property):
        lines = []
        for material in self.materials.keys():
            if property == 'material':
                for i in self.materials[material].keys():
                    lines.append(material)
            elif property == 'i':
                for i in self.materials[material].keys():
                    lines.append(i)
            else:
                for i in self.materials[material].keys():
                    lines.append(self.materials[material][i][property])
        return lines

    def add_property(self, run, property, value, i, fxn=None):
        material = self.get_material(run)
        self.materials[material][i][property] = value
        if (property, fxn) not in self.properties:
            self.properties.append((property, fxn))

    def add_ts_property(self, run, property, value, i, fxn=None):
        material = self.get_material(run)
        self.materials[material][i][property] = value
        if (property, fxn) not in self.ts_properties:
            self.ts_properties.append((property, fxn))

    def add_material(self, base, start, final, ovac, ovac_final, ub, lb, i):
        material = self.get_material(base, start, final)
        if material not in self.materials:
            self.materials[material] = {}
        self.materials[material][i] = {}
        new_material = self.materials[material][i]
        new_material['ub'] = ub
        new_material['lb'] = lb
        new_material['ovac'] = ovac
        new_material['ovac_final'] = ovac_final
        for property, fxn in self.properties:
            if fxn:
                value = fxn(self, base, start, final)
                self.add_property(base, property, value, i, fxn=fxn)




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

#%% Plot
def get_band_center(vasprun, orbital, element='O'):
    (title, columns, scaling_factors) = Make_Dos.make_dos(vasprun, ['{}:{}'.format(orbital, element)])
    maximum = np.argmax(np.array(columns[0]) > 0)
    energies = columns[0]
    dos = np.array(columns[3]) - np.array(columns[4])
    return np.average(energies[:maximum], weights=dos[:maximum])

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

def get_bars(mat, ts, suffix='', style=''):
    lb_list = ts['lb{}'.format(suffix)]
    ub_list = ts['ub{}'.format(suffix)]

    if 'paths' not in mat:
        bulk_e = mat['base']['energy']
        # type(mat['start{}'.format(suffix)])
        start_bars = [mat['start{}'.format(suffix)]['energy'] if mat['start{}'.format(suffix)].__class__ == dict else 0]
        final_bars = [x['energy'] if x.__class__ == dict else 0 for x in mat['finals{}'.format(suffix)]]
        bars = start_bars + final_bars

        lb_bars = [0]
        delta = np.array([0])
        for i, energy in enumerate(bars):
            if i == 0:
                start_energy = energy
            else:
                lb_bars.append(
                    lb_list[i - 1]['energy'] - min(energy, start_energy) if energy < 0 and start_energy < 0 and
                        lb_list[i - 1].__class__ == dict else 0)
                delta = np.append(delta, np.array(abs(start_energy - energy)))

        ub_bars = [0]
        for i, energy in enumerate(bars):
            if i == 0:
                start_energy = energy
            else:
                if len(ub_list[i - 1]) >= 2:
                    ub_bars.append(max([max(x['energies']) for x in ub_list[i - 1]]) - min(energy,
                                                                                           start_energy) if energy < 0 and start_energy < 0 and len(
                        ub_list[i - 1]) > 0 and ub_list[i - 1][0].__class__ == dict else 0)
                else:
                    ub_bars.append(0)

        bars = [(x + get_FERE_chemical_potential('O')) - bulk_e if x != 0 else 0 for x in bars]
        return bars, lb_bars, ub_bars, delta
    else: # Get Thermodynamic information
        paths = mat['paths']
        bulk_e = mat['base']['energy']
        bars = []
        lb_bars = []
        ub_bars = []

        # Defining functions to get energy
        def get_energy(i):
            key = i
            if key in mat and mat[key]:
                return mat[key]['energy']
            else:
                return 0
        def get_ts_energies(path, energies):
            if path in ts['lb'] and ts['lb'][path]:
                lb = ts['lb'][path]['energy'] - min(energies)
            else:
                lb = 0
            if path in ts['ub'] and ts['ub'][path] and len(ts['ub'][path]) >= 2 and type(ts['ub'][path][0]) == dict:
                ub = max( [max(x['energies']) - min(energies) for x in ts['ub'][path]  ])
            else:
                ub = 0
            return lb, ub

        # get information for each path
        for path in paths:
            energies = (get_energy(path[0]), get_energy(path[1]))
            lb, ub = get_ts_energies(path, energies)
            bars.append(energies)
            lb_bars.append(lb)
            ub_bars.append(ub)

        return bars, lb_bars, ub_bars

def make_bar_chart(mat, bars_list, delta, width=0.75, colors=['red', 'blue', 'green', 'yellow'], figsize=(10,3), text=True):
    # xlabels = ['Start'] + ['Final.{}'.format(i) for i in range(len(bars) - 1)]
    xlabels = ['Start'] + [
        '{} through {}'.format(x['pathway_id']['initial-final'], x['pathway_id']['critical-triangle']) if type(
            x) == dict and 'pathway_id' in x else 'Final.{}'.format(i) for i, x in enumerate(mat['finals'])]
    legend = ['Ovac', 'LB Diffusion', 'UB Diffusion']
    width = .75 / len(bars_list)
    x = np.arange(len(bars_list[0]))
    fig, ax = plt.subplots(figsize=figsize)# type: plt.Figure, plt.Axes
    for i, bars in enumerate(bars_list):
        if i >= 1:
            # lower = np.array(bars) - delta
            lower = [bars[i] - delta[i] if bars[i] > 0 else 0 for i in range(len(bars))]
            delta_temp = [delta[i] if bars[i] > 0 else 0 for i in range(len(bars))]
            ax.bar(x + i * width, lower, width, color=colors[i])
            ax.bar(x + i * width, delta_temp, width, color='dark{}'.format(colors[i]), bottom=lower,
                           label='_nolegend_')
        else:
            ax.bar(x + i * width, bars, width, color=colors[i])

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
    ax.set_xticklabels(range(len(bars)))
    ax.set_ylabel('Energy (eV)')
    ax.set_title(mat['base']['material'])
    ax.legend(legend, bbox_to_anchor=(0.5, -0.1), loc='upper center')
    ax.set_xticks(x + width * len(bars_list) / 2)
    ax.set_xticklabels(xlabels)
    if text:
        for i in range(len(bars_list[0])):
            text_x = x[i] + 1*width - width/2
            text_y = float(bars_list[1][i])
            ax.text(text_x,text_y, '{:.3}'.format(text_y), fontdict={'size' : 'xx-small'})

            text_x = x[i] + 2*width - width/2
            text_y = float(bars_list[2][i])
            ax.text(text_x,text_y, '{:.3}'.format(text_y), fontdict={'size' : 'xx-small'})

    fig.show()
    return fig, ax

def plot_full(mat_dict, material, fig, ax, bars_list,
              labels = [], label_all=False, string='', plot_higher_barriers=False, colors=['red', 'blue', 'green'],
              fit=False, linewidth=1.5, max_diff=5, csv_dict:Csv_dict=None):
    ts_estimates = [np.average([bars_list[1][i], bars_list[2][i]]) if bars_list[2][i] else 9999 for i in range(len(bars_list[0]))]
    lowest_ts_estimate = min(ts_estimates)
    lowest_ovac = 9999
    for i in range(len(bars_list[0])):
        lb = bars_list[1][i]
        ub = bars_list[2][i]
        ovac_e = min(bars_list[0][i], bars_list[0][0])
        if ovac_e and lb and ub:
            if 'spinel' in mat_dict[material]['base']['material']:
                color = colors[0]
            elif 'perovskite' in mat_dict[material]['base']['material']:
                color = colors[1]
            else:
                color = colors[2]
            zorder = 2
            if (ts_estimates[i] > lowest_ts_estimate) or ub - lb > max_diff:
                if not plot_higher_barriers:
                    continue
                zorder = 1
                if color == 'red':
                    color = '#F5B7B1'
                    # color = 'red'
                elif color == 'blue':
                    color = '#AED6F1'
                    # color = 'blue'
                else:
                    color = 'light' + color
            else:
                lowest_ovac = ovac_e
                if material in labels or label_all:
                    if ub - lb < max_diff:
                        ax.text(ovac_e, lb, material)
                if fit:
                    string, predicted = add_to_string(string, mat_dict[material]['base'], ovac_e,
                                                      (lb, ub),
                                                      weight=3, max_val=1.5)
                    error_min.append(min(abs(lb - predicted), abs(predicted - ub)))
                    error_max.append(max(abs(lb - predicted), abs(predicted - ub)))
                    ax_pred.plot([predicted, predicted], [lb, ub], marker='_',
                                 color=color,
                                 linewidth=linewidth, markersize=8, markeredgewidth=linewidth, zorder=zorder)
            #                     ax.plot([bars_list[0][i]], [ predicted ], marker='o', color=color, markeredgecolor=color, markerfacecolor='none',
            #                              linewidth=linewidth, markersize=5, markeredgewidth=linewidth, zorder=zorder)

            # if ub - lb < max_diff:
            if True:
                if csv_dict:
                    raise Exception('Dont add to csv_dict')
                    csv_dict.add_material(mat_dict[material]['base'], mat_dict[material]['start'],
                                          mat_dict[material]['finals'][i-1], ovac_e, ub, lb)
                ax.plot([ovac_e, ovac_e], [lb, ub], marker='_',
                        color=color,
                        linewidth=linewidth, markersize=8, markeredgewidth=linewidth, zorder=zorder)
    return fig, ax, lowest_ts_estimate, lowest_ovac

