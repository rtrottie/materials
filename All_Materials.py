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
from pymatgen.io.lobster import Icohplist
from Database_Tools import get_file
from Helpers import get_FERE_chemical_potential
from scipy.stats import gmean
import numpy as np
import Make_Dos
from pymatgen.io.ase import AseAtomsAdaptor
import time
import functools
from pymatgen.analysis.local_env import VoronoiNN

from Make_Dos import make_dos
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
    'ts_type' : {'$exists' : False},
    'material' : {'$nin' : ['h-barrier', 'from_zach', 'zunger']}
}

materials = [x['material'] for x in db.database.find(match_criteria).sort('material')]

if verbose:
    print(len(materials))
    print(materials)


def get_matches(temp_match_criteria, find_highest=False):
    temp_matches = list(db.database.find(temp_match_criteria).sort('energy'))
    key = 'defect_type'
    if len(temp_matches) > 1:
        best_match = temp_matches[0]
        if find_highest:
            nups = { run['incar']['NUPDOWN']: run for run in temp_matches } # initialize dictionary
            for run in temp_matches: # find highest energy run for each nupdown
                if nups[run['incar']['NUPDOWN']]['energy'] < run['energy']:
                    nups[run['incar']['NUPDOWN']] = run
            return nups[functools.reduce(lambda a,b : a if nups[a]['energy'] < nups[b]['energy'] else b, nups.keys())]# find lowest energy nup
            # for nup in nups:
            #     if match['incar']['NUPDOWN'] == best_match['incar']['NUPDOWN']:
            #         best_match = match

        #         view_multiple(temp_matches)
        # print('Multiple Matches {}'.format(temp_match_criteria))
        # [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        # [print('{} {}\n'.format(x[key], x['energy'])) for x in temp_matches if key in x]
        # print()
        return best_match
    elif len(temp_matches) == 0:
        # print('No Matches {}'.format(temp_match_criteria))
        return None
    else:
        return temp_matches[0]


mat_dict = {}
ts_dict = {}

for material in materials:
    # print(material)
    mat_dict[material[1]] = {}
    ts_dict[material[1]] = {}
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
    match_criteria = {
        'defect_location': {'$exists': False},
        'material': material,
        'labels': {'$nin': ['surface', 'ts', 'interpolation', 'charged_defect']},
        'pathway_count': {'$exists': True},
        'poscar.structure.lattice.a': {'$lt': l},
        'ts_type': {'$exists': False},
        'energy' : {'$exists' : True},
        'misc_labels' : {'$nin' : ['kpts_div_2']},
    }
    unit = get_matches(match_criteria)
    lb_criteria = {
        'material': material,
        'ts_type': {'$all' : ['pc']},
        'index': -1,
        'misc_labels' : {'$nin' : ['kpts_div_2']},
        'energy' : {'$exists' : True},
        'kpoints.kpoints' : base['kpoints']['kpoints'],
        'labels': {'$nin': ['surface', 'interpolation', 'charged_defect']},
    }
    ub_criteria = {
        'material': material,
        'interpolation_direction': {'$exists': True},
        'index': -1,
        'misc_labels' : {'$nin' : ['kpts_div_2']},
        'kpoints.kpoints' : base['kpoints']['kpoints'],
        'labels': {'$nin': ['surface', 'charged_defect']},
    }

    match_criteria = {
        'defect_location': 'start',
        'material': material,
        'misc_labels' : {'$nin' : ['kpts_div_2']},
        'energy' : {'$exists' : True},
        'kpoints.kpoints' : base['kpoints']['kpoints'],
        'labels': {'$nin': ['surface', 'ts', 'interpolation', 'charged_defect']},
    }
    start = get_matches(match_criteria)

    mat_dict[material[1]]['unit'] = unit
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

        lb = get_matches(lb_criteria, True)
        ub = list(db.database.find(ub_criteria))
        interp_directions = [run['interpolation_direction'][0] for run in ub]
        # print(interp_directions)
        continue_checking = True
        checks = ['1', '2']
        while continue_checking:
            continue_checking = False
            for check in checks.copy():
                if check + '.1' and check + '.2' in interp_directions:
                    continue_checking = True
                    checks.remove(check)
                    checks.append(check + '.1')
                    checks.append(check + '.2')
        # print('check')
        # print(interp_directions)
        # print(checks)
        cleaned_ubs = {}
        for run in ub:
            id = run['interpolation_direction'][0]
            if id not in checks:
                continue
            if id not in cleaned_ubs:
                cleaned_ubs[id] = run
            else:
                if max(run['energies']) < max(cleaned_ubs[id]['energies']):
                    cleaned_ubs[id] = run
        ub = [cleaned_ubs[id] for id in cleaned_ubs.keys()]

        mat_dict[material[1]]['finals'].append(final)
        ts_dict[material[1]]['lb'].append(lb)
        ts_dict[material[1]]['ub'].append(ub)




#%% Plotting Functions

string = ''

fit = True

done = [
    'bacoo3',
    'bafeo3',
    'bahfo3',
    'bamn3po3',
    'bicoo3',
    'biino3',
    'bisbo3',

    'basno3',
    'baseo3',
    'batao3',
    'bazro3',
    'bialo3',
    'bimno3',
    'bisco3',
    'bivo3',
    'cahfo3',
    'camno3',
    'catio3',
    'cacro3',
    'caseo3',
    # 'casno3',
    'cavo3',
    'inalo3',
    'kbio3',
    'ksbo3',
    'ktao3',
    'lacro3',
    'laalo3',
    'lagao3',
    'laino3',
    'lamn3po3',
    'lasco3',
    'latio3',
    'mgseo3',
    # 'mgcr2o4',
    'natao3',
    'sncoo3',
    'snnio3',
    'sntao3',
    'sralo3',
    'srcoo3',
    'srfeo3',
    'srhfo3',
    'srmn3po3',
    'srmno3',
    'srseo3',
    'srsno3',
    'srtio3',
    # 'ygao3',
    'ycro3',

    'nabio3',
    'snhfo3',
    'srvo3',
    'yalo3',
    'ygao3',
    'ynio3',
    'scalo3',
    'snmno3',
    'srseo3',
    'srzro3',
    'mggeo3',
    'mgmno3',
    'zncro3',
    'zngeo3',
    'znmno3',
    'nasbo3',
    #
    'cdal2o4',
    'coal2o4',
    'coga2o4',
    'cral2o4',
    'cual2o4',
    # 'feal2o4',
    'femn2o4',
    'gemg2o4',
    'geni2o4',
    'mgal2o4',
    'mgcr2o4',
    'mnal2o4',
    # 'nial2o4',
    'vmg2o4',
    'nial2o4',
    'znal2o4',

]

done_lower_ub = [
    'biino3',
    'casno3',
    'geni2o4',
    'lasco3',
    'mgseo3',
    'sncoo3',
    'snnio3',
    'zngeo3',
    'laino3',
    'mnseo3',
    'zncro3',
    'coseo3',
    'lacro3',
    'srseo3',
    'snseo3',
    'snzro3',
]

converging = [
    'lavo3',
    'mgcro3',
    'zncoo3',
    'znseo3',
    'lamno3',
    'mgsio3',
    'cucr2o4',
    'cafeo3',
    'fesio3',
    'layo3',
    'limn2o4',
    'mnsio3',

]

bad_energies = [
    'vsio3', # bad vacancies
    'bicuo3', # bad vacancies
    'bifeo3', # bad vacancies
    'binio3', # bad vacancies
    'crcr2o4', # bad vacancies
    'mnv2o4', # bad vacancies
    'tisio3', # bad vacancies
    'snbeo3', # bad vacancies
    'tibeo3', # bad vacancies
    'zrbeo3', # bad vacancies
    'bigao3', # Bad LB NEB
    'nife2o4', # BAD LB
    'bicro3', # Bad vacancies and LB
    'yfeo3', # Bad vacancies and LB
    'gefe2o4', # BAD UB
    'gesio3', #BAD UB
    'nini2o4', # BAD UB
    'snsio3', #BAD UB
    'alv2o4', # NO SPIN
    'navo3', # Deletetd
    'yvo3', # deleted
]

bad_ub = [
    'lacro3',
]

nup_wrong = [
]

lookatit = [

]


checked = [
]

try2 = [
]

failed_making_neb = [
    'snseo3',
    'snzro3',
    'alal2o4',
    # 'alv2o4',
    'canio3',
    'coseo3', # bad energy
    'fecr2o4',
    'mnseo3',
    'snfeo3',
    'sncro3',
    'sntio3',
    'tiseo3',
    'ytao3',
    'sicd2o4',
    'cuga2o4',
    'snvo3',


]

misaligned = [
]

no_nupdown = [
]
gone = [
    'vfe2o4',
    'tife2o4',
    'tico2o4',
    'sizn2o4',
    'simg2o4',
    'cocr2o4',
    'crcr2o4'
    'femn2o4',
    'mnmn2o4',

]
plot = [
    'BaSnO3',
    'BaHfO3',
    'BaZrO3',
    'CaSnO3',
    'CaHfO3',
    'YCrO3',
    'LaCrO3',
    'LaScO3',
    'LaAlO3',
    'NaTaO3',
    'KTaO3',
    'CaTiO3',
    'SrTiO3',
    'SrZrO3',
    'SrHfO3',
    'SrSnO3',
]

plot = [x.lower() for x in plot]

aluminates = [
]
no_vacanciestxt = [
]
#%% Properties



def get_fit_info(base, pathways, lbs, ubs):
    for i in range(len(pathways)):
        start, final = pathways[i]
        lb = lbs[i]
        ub = ubs[i]
        if not (start and final and lb and ub): # check if everything exists
            continue


def find(site, sites, info=False):
    for i in range(len(sites)):
        if info:
            print(sites[i].distance(site))
        if sites[i].distance(site) < 1e-8:
            return i
    raise Exception('Could Not Find')


def get_coords(perfect, defect):
    s_perfect = perfect.copy()
    s_defect = defect.copy()
    for site in s_defect:
        s_perfect.append(site.specie, site.frac_coords)
    s = s_perfect.copy() # Combine structures
    dist = 0.25

    while len(s) > 2:
        site_is = [ find(site[0], s) for site in s.get_sites_in_sphere(s[0].coords, dist) ]
        if len(site_is) == 2 and s[site_is[0]].specie == s[site_is[1]].specie:
            s.remove_sites(site_is)
        try:
            site_is = [ find(site[0], s) for site in s.get_sites_in_sphere(s[1].coords, dist) ]
            if len(site_is) == 2 and s[site_is[0]].specie == s[site_is[1]].specie:
                s.remove_sites(site_is)
        except:
            dist += 0.001
            pass
        dist += 0.001
        # print(len(s))
    if len(s) > 1 and s[0].distance(s[1]) > 1:
        raise Exception('Site is not the same')
    return find(s[0], perfect)


def get_vasprun(run):
    v = get_file(fs, run['vasprun'], 'vasprun')
    vasprun = Vasprun(v)
    os.remove(v)
    return vasprun


class FitProperty:
    def __init__(self, name, run='unit'):
        self.name = name
        self.run = run

    def get_run(self, unit, base, start, final):
        if self.run == 'unit':
            run = unit
        elif self.run == 'base':
            run = base
        elif self.run == 'start':
            run = start
        elif self.run == 'final':
            run = final
        else:
            raise Exception('Invalid run provided: {}'.format(self.run))
        return run

    def get_from_run(self, run):
        if self.name in run:
            return run[self.name]
        else:
            return None

    def fxn(self, unit, base, start, final):
        Exception('Need To Define')

    do_update = False
    force_update = False
    print_update = True
    hard_force_update = False

    def set_update(self):
        self.do_update = True
        self.force_update = True
        self.print_update = False

    def update(self, run, property):
        if self.print_update:
            print('{} : {} {}'.format(run['material'], self.name, property))
        if self.do_update:
            run[self.name] = property
            db.database.update_one({'_id': run['_id']}, {'$set': {self.name: property}})
            pass

    def get(self, unit, base, start, final):
        run = self.get_run(unit, base, start, final)
        if not self.force_update and not self.hard_force_update:
            property = self.get_from_run(run)
            if property != None:
                return property
        property = self.fxn(unit, base, start, final)
        self.update(run, property)
        return property


class TSProperty(FitProperty):

    def get_run(self, unit, base, start, final, lb, ub):
        if self.run == 'unit':
            run = unit
        elif self.run == 'base':
            run = base
        elif self.run == 'start':
            run = start
        elif self.run == 'final':
            run = final
        elif self.run == 'ub':
            run = ub
        elif self.run == 'lb':
            run = lb
        else:
            raise Exception('Invalid run provided: {}'.format(self.run))
        return run


    def get(self, unit, base, start, final, lb, ub):
        run = self.get_run(unit, base, start, final, lb, ub)
        if not self.force_update:
            property = self.get_from_run(run)
            if property != None:
                return property
        property = self.fxn(unit, base, start, final, lb, ub)
        self.update(run, property)
        return property

    def fxn(self, unit, base, start, final, lb, ub):
        Exception('Need To Define')


class BandCenter(FitProperty):
    do_update = True
    def __init__(self, name, bandcenter, run='unit'):
        self.name = name
        self.run = run
        self.bandcenter = bandcenter

    def fxn(self, unit, base, start, final):
        run = self.get_run(unit, base, start, final)
        vasprun = get_vasprun(run)
        labels, data, scale = make_dos(vasprun, self.bandcenter)
        energies = np.array(data[0])
        dos = np.array(data[-2]) - np.array(data[-1])
        max = np.where(energies<0)[0][-1] + 1
        bandcenter = np.average(energies[:max], weights=dos[:max])
        print(bandcenter)
        return bandcenter


class BandCenterWithConductionBand(BandCenter):
    do_update = True

    def fxn(self, unit, base, start, final):
        run = self.get_run(unit, base, start, final)
        vasprun = get_vasprun(run)
        labels, data, scale = make_dos(vasprun, self.bandcenter)
        energies = np.array(data[0])
        dos = np.array(data[-2]) - np.array(data[-1])
        bandcenter = np.average(energies, weights=dos)
        print(bandcenter)
        return bandcenter


class VoronoiInfo(FitProperty):
    force_update = False
    do_update = True

    def fxn(self, unit, base, start, final):
        material = get_lower_material(unit)
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))

        base_s = Poscar.from_dict(base['poscar']).structure
        start_s = Poscar.from_dict(start['poscar']).structure
        run_s = Poscar.from_dict(self.get_run(unit, base, start, final)['poscar']).structure
        # print(len(base_s))
        # print(len(start_s))
        start_i = get_coords(base_s.copy(), run_s.copy())
        print(start_i)
        vnn = VoronoiNN(targets=[Element(x) for x in elements])
        start_vnn = vnn.get_nn_info(base_s, start_i)

        weight = {}
        total_length = {}
        for pt in start_vnn:
            # print('{}: {:3.2f}'.format(pt['site_index'], pt['weight']))
            temp_weight = round(pt['weight']+0.45)
            if pt['site'].species_string in weight:
                weight[pt['site'].species_string] = weight[pt['site'].species_string] + temp_weight
            else:
                weight[pt['site'].species_string] = temp_weight
            if temp_weight >= 0.99:
                bond_length = base_s.get_distance(start_i, pt['site_index'])
                if pt['site'].species_string in total_length:
                    total_length[pt['site'].species_string] = total_length[pt['site'].species_string] + bond_length
                else:
                    total_length[pt['site'].species_string] = bond_length

        if material.find(elements[0].lower()) == 0 and material.find(elements[1].lower()) == 0:
            raise Exception('Both Could be Element A')
        elif material.find(elements[0].lower()) == 0:
            bonds = {
                'A_Bonds' : weight[elements[0]],
                'A_length' : total_length[elements[0]] / weight[elements[0]],
                'B_Bonds' : weight[elements[1]],
                'B_length' : total_length[elements[1]] / weight[elements[1]],
                'avg_length' : (total_length[elements[0]] + total_length[elements[1]]) / (weight[elements[0]] + weight[elements[1]]),
            }
        elif material.find(elements[1].lower()) == 0:
            bonds = {
                'A_Bonds' : weight[elements[1]],
                'A_length' : total_length[elements[1]] / weight[elements[1]],
                'B_Bonds' : weight[elements[0]],
                'B_length' : total_length[elements[0]] / weight[elements[0]],
                'avg_length' : (total_length[elements[0]] + total_length[elements[1]]) / (weight[elements[0]] + weight[elements[1]]),
            }
        else:
            raise Exception('Couldn\'t find A')
        return bonds


class Efermi(FitProperty):
    do_update = True
    def fxn(self, unit, base, start, final):
        run = self.get_run(unit, base, start, final)
        vasprun = get_vasprun(run)
        return vasprun.efermi


class Bandgap(FitProperty):
    do_update = True
    def fxn(self, unit, base, start, final):
        run = self.get_run(unit, base, start, final)
        vasprun = get_vasprun(run)
        bg = vasprun.get_band_structure().get_band_gap()['energy']
        return bg


class Volume(FitProperty):
    do_update = True
    def fxn(self, unit, base, start, final):
        s = Poscar.from_dict(unit['poscar']).structure
        return s.volume / s.num_sites


class BondStrength(FitProperty):
    do_update = True
    force_update = False

    def fxn(self, unit, base, start, final):
        f = get_file(fs, unit['ICOHPLIST_lobster'])
        icohp = Icohplist(filename=f).icohpcollection
        os.remove(f)
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))
        elements.append('O')  # Put O at end for easier referencing
        material = get_lower_material(unit)
        unit_s = Poscar.from_dict(unit['poscar']).structure

        summed_icohp = {}
        for i, site in enumerate(unit_s):
            icohp_dict = icohp.get_icohp_dict_of_site(i)
            for icohp_key in icohp_dict:
                icohp_value = icohp_dict[icohp_key].summed_icohp
                # print('{}: {:3.2f}'.format(pt['site_index'], pt['weight']))
                if site.species_string in summed_icohp:
                    summed_icohp[site.species_string] = (
                        summed_icohp[site.species_string][0] + icohp_value,
                        summed_icohp[site.species_string][1] + 1
                    )
                else:
                    summed_icohp[site.species_string] = (icohp_value, 1)

        if material.find(elements[0].lower()) == 0 and material.find(elements[1].lower()) == 0:
            raise Exception('Both Could be Element A')
        elif material.find(elements[0].lower()) == 0:
            bonds = {
                'A_Bonds' : summed_icohp[elements[0]][0] / summed_icohp[elements[0]][1],
                'B_Bonds' : summed_icohp[elements[1]][0] / summed_icohp[elements[1]][1]
            }
        elif material.find(elements[1].lower()) == 0:
            bonds = {
                'A_Bonds' : summed_icohp[elements[1]][0] / summed_icohp[elements[1]][1],
                'B_Bonds' : summed_icohp[elements[0]][0] / summed_icohp[elements[0]][1]
            }
        else:
            raise Exception('Couldn\'t find A')
        bonds['O_bonds'] = summed_icohp['O'][0] / summed_icohp['O'][1]
        return bonds


class DiffusionInfo(FitProperty):


    do_update = True
    def get_position_along_path(self, start_coords, final_coords, dist):
        '''
        Interpolate between start and final
        :param start_coords:
        :param final_coords:
        :param dist: between 0 and 1
        :return:
        '''
        start_coords = np.array(start_coords)
        final_coords = np.array(final_coords)
        vector = final_coords - start_coords
        [x + 1 if x < -0.5 else x - 1 if x > 0.5 else x for x in vector]  # Make vector between -0.5 and 0.5
        vector = vector * dist
        return start_coords + vector

    def get_nearest_atoms(self, base_s: Structure, start_s: Structure, final_s: Structure):
        '''
        Gets IDPP interpolated pathway between start and final
        :param base_s: Structure
        :param start_s: Structure
        :param final_s: Structure
        :return: list(Structure)
        '''

        start_i = get_coords(base_s.copy(), start_s.copy())
        final_i = get_coords(base_s.copy(), final_s.copy())

        start_atom = base_s[start_i]  # type: PeriodicSite
        final_atom = base_s[final_i]  # type: PeriodicSite

        center_location = self.get_position_along_path(start_atom.frac_coords, final_atom.frac_coords, 0.5)

        close_atoms = base_s.get_sites_in_sphere(base_s.lattice.get_cartesian_coords(center_location), 5, include_image=True, include_index=True)
        close_atoms.sort(key=lambda x: x[1])
        close_atoms = [ x for x in close_atoms if x[2] not in [start_i, final_i] and x[0].specie != Element('O')]
        close_atoms = close_atoms[:5]
        return close_atoms, start_i, final_i, center_location


    def get_sorted_nearest_bonds(self, structure):
        return None


    def fxn(self, unit, base, start, final):
        '''
        Get Area of Critical Triangle
        :param unit:
        :param base:
        :param start:
        :param final:
        :return:
        '''

        base_s = Poscar.from_dict(base['poscar']).structure
        start_s = Poscar.from_dict(start['poscar']).structure
        final_s = Poscar.from_dict(final['poscar']).structure
        close_atoms, start_i, final_i, center_location = self.get_nearest_atoms(base_s, start_s, final_s)
        coords = [x[0].coords for x in close_atoms]

        A = np.sum(np.square(coords[1] - coords[0]))
        B = np.sum(np.square(coords[2] - coords[0]))
        C = np.sum(np.square(coords[1] - coords[2]))

        S = np.sqrt(4*A*B - (C - A - B)**2)/4

        vis = base_s.copy()
        vis.replace(start_i, Element('H'))
        vis.replace(final_i, Element('H'))
        vis.append(Element('H'), center_location)

        # view_multiple([AseAtomsAdaptor.get_atoms(vis)], is_ase=True)
        # time.sleep(0.5)
        return S


class DiffusionAdjustedInfo(DiffusionInfo):

    do_update = True

    def fxn(self, unit, base, start, final):
        '''
        Get Area of Critical Triangle
        :param unit:
        :param base:
        :param start:
        :param final:
        :return:
        '''

        base_s = Poscar.from_dict(base['poscar']).structure
        start_s = Poscar.from_dict(start['poscar']).structure
        final_s = Poscar.from_dict(final['poscar']).structure
        close_atoms, start_i, final_i, center_location = self.get_nearest_atoms(base_s, start_s, final_s)
        coords = [x[0].coords for x in close_atoms]

        material = get_lower_material(unit)
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))
        if material.find(elements[0].lower()) == 0:
            A_species = elements[0]
            B_species = elements[1]
        else:
            A_species = elements[1]
            B_species = elements[0]

        al = np.arccos(np.dot(coords[2]-coords[0], coords[1]-coords[0]) / np.linalg.norm(coords[2]-coords[0]) / np.linalg.norm(coords[1]-coords[0]))
        be = np.arccos(np.dot(coords[2]-coords[1], coords[0]-coords[1]) / np.linalg.norm(coords[2]-coords[1]) / np.linalg.norm(coords[0]-coords[1]))
        ga = np.arccos(np.dot(coords[0]-coords[2], coords[1]-coords[2]) / np.linalg.norm(coords[0]-coords[2]) / np.linalg.norm(coords[1]-coords[2]))

        bonded_area = 0
        for i, angle in enumerate([al, be, ga]):
            if close_atoms[i][0].species_string == A_species:
                bonded_area += (final['bond']['A_length']/2)**2 * angle / np.pi
            elif close_atoms[i][0].species_string == B_species:
                bonded_area += (final['bond']['B_length']/2)**2 * angle / np.pi
            else:
                raise Exception('Elements inccorecttly identified')


        A = np.sum(np.square(coords[1] - coords[0]))
        B = np.sum(np.square(coords[2] - coords[0]))
        C = np.sum(np.square(coords[1] - coords[2]))

        S = np.sqrt(4*A*B - (C - A - B)**2)/4

        # vis = base_s.copy()
        # vis.replace(start_i, Element('H'))
        # vis.replace(final_i, Element('H'))
        # vis.append(Element('H'), center_location)

        # view_multiple([AseAtomsAdaptor.get_atoms(vis)], is_ase=True)
        # time.sleep(0.5)
        return S - bonded_area


class DiffusionAdjustedFullInfo(DiffusionInfo):
    force_update = False
    do_update = True

    def fxn(self, unit, base, start, final):
        '''
        Get Area of Critical Triangle
        :param unit:
        :param base:
        :param start:
        :param final:
        :return:
        '''

        base_s = Poscar.from_dict(base['poscar']).structure
        start_s = Poscar.from_dict(start['poscar']).structure
        final_s = Poscar.from_dict(final['poscar']).structure
        close_atoms, start_i, final_i, center_location = self.get_nearest_atoms(base_s, start_s, final_s)
        coords = [x[0].coords for x in close_atoms]

        material = get_lower_material(unit) #type: str
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))
        if material.find(elements[0].lower()) == 0:
            A_species = elements[0]
            B_species = elements[1]
        else:
            A_species = elements[1]
            B_species = elements[0]

        al = np.arccos(np.dot(coords[2]-coords[0], coords[1]-coords[0]) / np.linalg.norm(coords[2]-coords[0]) / np.linalg.norm(coords[1]-coords[0]))
        be = np.arccos(np.dot(coords[2]-coords[1], coords[0]-coords[1]) / np.linalg.norm(coords[2]-coords[1]) / np.linalg.norm(coords[0]-coords[1]))
        ga = np.arccos(np.dot(coords[0]-coords[2], coords[1]-coords[2]) / np.linalg.norm(coords[0]-coords[2]) / np.linalg.norm(coords[1]-coords[2]))

        bonded_area = 0
        for i, angle in enumerate([al, be, ga]):
            if close_atoms[i][0].species_string == A_species:
                bonded_area += (final['bond']['A_length'])**2 * angle / np.pi
            elif close_atoms[i][0].species_string == B_species:
                bonded_area += (final['bond']['B_length'])**2 * angle / np.pi
            else:
                raise Exception('Elements inccorecttly identified')


        A = np.sum(np.square(coords[1] - coords[0]))
        B = np.sum(np.square(coords[2] - coords[0]))
        C = np.sum(np.square(coords[1] - coords[2]))

        S = np.sqrt(4*A*B - (C - A - B)**2)/4

        # vis = base_s.copy()
        # vis.replace(start_i, Element('H'))
        # vis.replace(final_i, Element('H'))
        # vis.append(Element('H'), center_location)

        # view_multiple([AseAtomsAdaptor.get_atoms(vis)], is_ase=True)
        # time.sleep(0.5)
        return S - bonded_area


class DiffusionDistance(FitProperty):
    force_update = False
    do_update = True

    def fxn(self, unit, base, start, final):
        '''
        Get Area of Critical Triangle
        :param unit:
        :param base:
        :param start:
        :param final:
        :return:
        '''

        base_s = Poscar.from_dict(base['poscar']).structure
        start_s = Poscar.from_dict(start['poscar']).structure
        final_s = Poscar.from_dict(final['poscar']).structure

        start_i = get_coords(base_s.copy(), start_s.copy())
        final_i = get_coords(base_s.copy(), final_s.copy())

        distance = base_s.get_distance(start_i, final_i)
        return distance


class FormationEnergy(FitProperty):
    do_update = True
    force_update = True
    chem_pots = {
        'Se' : -3.50402579890625,
    }

    def fxn(self, unit, base, start, final):
        unit_s = Poscar.from_dict(unit['poscar']).structure
        elemental_energies = 0
        for atom in unit_s:
            element = atom.species_string
            if element not in self.chem_pots:
                chem_pot_run = list(db.database.find({'material': {'$all': [element.lower(), 'watersplitting-chemical-potential']}}).sort('energy'))[0]
                chem_pot_s = Poscar.from_dict(chem_pot_run['poscar']).structure
                self.chem_pots[element] = chem_pot_run['energy'] / len(chem_pot_s)
            elemental_energies += self.chem_pots[element]
        # print(self.chem_pots)
        return (unit['energy'] - elemental_energies) / len(unit_s)


class AtomizationEnergy(FitProperty):
    do_update = True
    force_update = True
    print_update = False
    chem_pots = {
    }

    def fxn(self, unit, base, start, final):
        unit_s = Poscar.from_dict(unit['poscar']).structure
        elemental_energies = 0
        for atom in unit_s:
            element = atom.species_string
            if element not in self.chem_pots:
                chem_pot_run = list(db.database.find({'material': {'$all': [element.lower(), 'watersplitting-atomization-energy']}}).sort('energy'))[0]
                chem_pot_s = Poscar.from_dict(chem_pot_run['poscar']).structure
                self.chem_pots[element] = chem_pot_run['energy'] / len(chem_pot_s)
            elemental_energies += self.chem_pots[element]
        return -(unit['energy'] - elemental_energies) / len(unit_s)


class PymgDifference(DiffusionInfo):
    do_update = True
    force_update = False
    def fxn(self, unit, base, start, final):
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))

        base_s = Poscar.from_dict(base['poscar']).structure
        run_s = Poscar.from_dict(self.get_run(unit, base, start, final)['poscar']).structure
        start_i = get_coords(base_s.copy(), run_s.copy())
        print(start_i)
        vnn = VoronoiNN(targets=[Element(x) for x in elements])
        start_vnn = vnn.get_nn_info(base_s, start_i)

        pauling_total = 0
        count = 0
        for pt in start_vnn:
            # print('{}: {:3.2f}'.format(pt['site_index'], pt['weight']))
            if pt['weight'] >= 0.05:
                pauling_total += self.pauling[pt['site'].species_string]
                count += 1
        pauling_diff = self.pauling['O'] - pauling_total / count
        return pauling_diff


class PaulingDifference(DiffusionInfo):
    do_update = True
    force_update = False
    print_update = False
    pauling = {'H':2.2,'Li':0.98,'Be':1.57,'B':2.04,'C':2.55,'N':3.04,'O':3.44,'F':3.98,'Na':0.93,'Mg':1.31,'Al':1.61,
               'Si':1.9,'P':2.19,'S':2.58,'Cl':3.16,'K':0.82,'Ca':1,'Sc':1.36,'Ti':1.54,'V':1.63,'Cr':1.66,'Mn':1.55,
               'Fe':1.83,'Co':1.88,'Ni':1.91,'Cu':1.9,'Zn':1.65,'Ga':1.81,'Ge':2.01,'As':2.18,'Se':2.55,'Br':2.96,
               'Kr':3,'Rb':0.82,'Sr':0.95,'Y':1.22,'Zr':1.33,'Nb':1.6,'Mo':2.16,'Tc':1.9,'Ru':2.2,'Rh':2.28,'Pd':2.2,
               'Ag':1.93,'Cd':1.69,'In':1.78,'Sn':1.96,'Sb':2.05,'Te':2.1,'I':2.66,'Xe':2.6,'Cs':0.79,'Ba':0.89,
               'La':1.1,'Ce':1.12,'Pr':1.13,'Nd':1.14,'Sm':1.17, 'Eu': 1.2,'Gd':1.2,'Dy':1.22,'Ho':1.23,'Er':1.24,'Tm':1.25,
               'Lu':1.27,'Hf':1.3,'Ta':1.5,'W':2.36,'Re':1.9,'Os':2.2,'Ir':2.2,'Pt':2.28,'Au':2.54,'Hg':2,'Tl':1.62,
               'Pb':2.33,'Bi':2.02,'Po':2,'At':2.2,'Ra':0.9,'Ac':1.1,'Th':1.3,'Pa':1.5,'U':1.38,'Np':1.36,'Pu':1.28,
               'Am':1.3,'Cm':1.3,'Bk':1.3,'Cf':1.3,'Es':1.3,'Fm':1.3,'Md':1.3,'No':1.3}


    def fxn(self, unit, base, start, final):
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))

        base_s = Poscar.from_dict(base['poscar']).structure
        run_s = Poscar.from_dict(self.get_run(unit, base, start, final)['poscar']).structure
        start_i = get_coords(base_s.copy(), run_s.copy())
        print(start_i)
        vnn = VoronoiNN(targets=[Element(x) for x in elements])
        start_vnn = vnn.get_nn_info(base_s, start_i)

        pauling_total = 0
        count = 0
        for pt in start_vnn:
            # print('{}: {:3.2f}'.format(pt['site_index'], pt['weight']))
            if pt['weight'] >= 0.05:
                pauling_total += self.pauling[pt['site'].species_string]
                count += 1
        pauling_diff = self.pauling['O'] - pauling_total / count
        return pauling_diff


class PaulingDifferenceWeighted(DiffusionInfo):
    do_update = True
    force_update = False
    print_update = False
    pauling = {'H':2.2,'Li':0.98,'Be':1.57,'B':2.04,'C':2.55,'N':3.04,'O':3.44,'F':3.98,'Na':0.93,'Mg':1.31,'Al':1.61,
               'Si':1.9,'P':2.19,'S':2.58,'Cl':3.16,'K':0.82,'Ca':1,'Sc':1.36,'Ti':1.54,'V':1.63,'Cr':1.66,'Mn':1.55,
               'Fe':1.83,'Co':1.88,'Ni':1.91,'Cu':1.9,'Zn':1.65,'Ga':1.81,'Ge':2.01,'As':2.18,'Se':2.55,'Br':2.96,
               'Kr':3,'Rb':0.82,'Sr':0.95,'Y':1.22,'Zr':1.33,'Nb':1.6,'Mo':2.16,'Tc':1.9,'Ru':2.2,'Rh':2.28,'Pd':2.2,
               'Ag':1.93,'Cd':1.69,'In':1.78,'Sn':1.96,'Sb':2.05,'Te':2.1,'I':2.66,'Xe':2.6,'Cs':0.79,'Ba':0.89,
               'La':1.1,'Ce':1.12,'Pr':1.13,'Nd':1.14,'Sm':1.17, 'Eu': 1.2,'Gd':1.2,'Dy':1.22,'Ho':1.23,'Er':1.24,'Tm':1.25,
               'Lu':1.27,'Hf':1.3,'Ta':1.5,'W':2.36,'Re':1.9,'Os':2.2,'Ir':2.2,'Pt':2.28,'Au':2.54,'Hg':2,'Tl':1.62,
               'Pb':2.33,'Bi':2.02,'Po':2,'At':2.2,'Ra':0.9,'Ac':1.1,'Th':1.3,'Pa':1.5,'U':1.38,'Np':1.36,'Pu':1.28,
               'Am':1.3,'Cm':1.3,'Bk':1.3,'Cf':1.3,'Es':1.3,'Fm':1.3,'Md':1.3,'No':1.3}


    def fxn(self, unit, base, start, final):
        elements = unit['elements'].copy()
        elements.pop(elements.index('O'))

        base_s = Poscar.from_dict(base['poscar']).structure
        run_s = Poscar.from_dict(self.get_run(unit, base, start, final)['poscar']).structure
        start_i = get_coords(base_s.copy(), run_s.copy())
        print(start_i)
        vnn = VoronoiNN(targets=[Element(x) for x in elements])
        start_vnn = vnn.get_nn_info(base_s, start_i)

        pauling_total = 0
        weight = 0
        for pt in start_vnn:
            # print('{}: {:3.2f}'.format(pt['site_index'], pt['weight']))
            pauling_total += self.pauling[pt['site'].species_string] * pt['weight']
            weight += pt['weight']
        pauling_diff = self.pauling['O'] - pauling_total / weight
        return pauling_diff

class TSEFermiShifts(TSProperty):
    do_update = True
    def fxn(self, unit, base, start, final, lb, ub):
        run = self.get_run(unit, base, start, final, lb, ub)
        ts_efermi = get_vasprun(run).efermi
        return ts_efermi



properties = [
    ('o_p_band_center', BandCenter('o_p_band_center', ['O:p'])),
    # ('o_p_band_center_with_conduction', BandCenterWithConductionBand('o_p_band_center_with_conduction', ['O:p'])),
    ('band_gap', Bandgap('band_gap')),
    # ('band_gap_base', Bandgap('band_gap_base', run='base')),
    ('efermi', Efermi('efermi')),
    # ('efermi_vac', Efermi('efermi_vac', run='start')),
    # ('volume_per_atom', Volume('volume_per_atom', run='unit')),
    # ('bond', VoronoiInfo('bond', run='start')),
    # ('bond_final', VoronoiInfo('bond', run='final')),
    # ('crit_triangle', DiffusionInfo('crit_triangle', run='final')),
    # ('diffusion_distance', DiffusionDistance('diffusion_distance', run='final')),
    # ('lobster_icohp', BondStrength('lobster_icohp', run='unit')),
    # ('crit_triangle_adjusted', DiffusionAdjustedInfo('crit_triangle_adjusted', run='final')),
    # ('crit_triangle_adjusted', DiffusionAdjustedFullInfo('crit_triangle_full_adjusted', run='final')),
    ('formation_energy', FormationEnergy('formation_energy', run='unit')),
    ('atomization_energy', AtomizationEnergy('atomization_energy', run='unit')),
    # ('pauling_diff', PaulingDifference('pauling_diff', run='final')),
    # ('pauling_diff_start', PaulingDifference('pauling_diff', run='start')),
    # ('pauling_diff_weighted', PaulingDifferenceWeighted('pauling_diff_weighted', run='final')),
    # ('pauling_diff_start_weighted', PaulingDifferenceWeighted('pauling_diff_weighted', run='start')),
]

ts_properties = [
    ('ts_efermi', TSEFermiShifts('ts_efermi', run='lb'))
]

# properties = []
# ts_properties = []
# plot = []


#%% Plot All


plot_all = False
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
fig_all, ax_all = plt.subplots(figsize=(2.5, 2.5))
fig_pred, ax_pred = plt.subplots(figsize=(6, 5))
ax_pred.set_ylim([0, 4])
ax_pred.set_xlim([0, 4])

ax_all.plot([-1], [-1], color='blue', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.plot([-1], [-1], color='red', marker='_', linewidth=linewidth, markersize=8, markeredgewidth=linewidth, )
ax_all.legend(['Perovskites', 'Spinels'], loc='upper left')
ax_all.set_ylabel('Diffusion Barrier (eV)')
ax_all.set_xlabel('Oxygen Vacancy Energy (eV)')

ax_all.set_ylim([0,5])
ax_all.set_xlim([-2, 8])

dont_plot = done + done_lower_ub + checked + try2 + misaligned + bad_energies + converging + nup_wrong + no_nupdown + no_vacanciestxt + gone
dont_plot = done + checked + try2 + misaligned + bad_energies + converging + nup_wrong + no_nupdown + no_vacanciestxt + gone + failed_making_neb

labels = []
unit_criteria = {
            'defect_type': {'$exists': False},
            'ts_type': {'$exists': False},
            'material': {'$all': [], '$in': ['spinel', 'perovskite'], '$nin': ['from_zach']},
            'labels' : {
                '$all' : ['unit_cell'],
                '$nin' : ['surface']}
        }

count = 0
csv_dict = Csv_dict()
if __name__ == '__main__':
# if True:
    for material in list(mat_dict.keys())[:]:
      if material not in ['ceo2']:
      # if material in done:
      # if material in ['bialo3']
    # for material in ['bacoo3', 'bafeo3']:

        mat_unit_criteria = copy.deepcopy(unit_criteria)
        mat_unit_criteria['material']['$all'].append(material)
        units = list(db.database.find(mat_unit_criteria).sort('energy'))
        unit = units[0]
        mat = mat_dict[material]

        # vasprun = get_vasprun(unit)
        thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
        bars_list = [thermo, lb_bars, ub_bars]
        plot_full(mat_dict, material, fig_all, ax_all, bars_list, labels=labels, label_all=label_all, string=string,
                  max_diff=0.75, linewidth=linewidth, plot_higher_barriers=True)
        if plot_all or material in plot:
            make_bar_chart(mat, bars_list, delta)
        elif not plot_all and material not in dont_plot:
            make_bar_chart(mat, bars_list, delta)
        bulk = mat['base']
        start = mat['start']
        base = mat['base']
        finals = mat['finals']
        count = count + len(finals)
        for i in range(1, len(mat['finals'])+1):
            if thermo[i] and lb_bars[i] and ub_bars[i]:
                ub = ts_dict[material]['ub'][i-1]
                ub_e = max(max([x['energies'] for x in ub])) - start['energy']
            elif thermo[i] and lb_bars[i]:
                ub_e = None
            else:
                continue

            final = finals[i-1]
            lb = ts_dict[material]['lb'][i-1]

            ovac_e = start['energy'] - bulk['energy']
            ovac_final = start['energy'] - bulk['energy']
            lb_e = lb['energy'] - start['energy']
            csv_dict.add_material(mat_dict[material]['base'], mat_dict[material]['start'],
                                  mat_dict[material]['finals'][i - 1], ovac_e, ovac_final, ub_e, lb_e, i)
            print(1)
            print(material)
            print(unit['material'])
            print(start['material'])
            for name, property in properties:
                value = property.get(unit, base, start, finals[i-1])
                if type(value) == dict:
                    for k in value:
                        csv_dict.add_property(unit, '{}_{}'.format(name,k), value[k], i)
                else:
                    csv_dict.add_property(unit, name, value, i)

            for name, property in ts_properties:
                value = property.get(unit, base, start, finals[i-1], lb, ub)
                if type(value) == dict:
                    for k in value:
                        csv_dict.add_property(unit, '{}_{}'.format(name,k), value[k], i)
                else:
                    csv_dict.add_ts_property(unit, name, value, i)
    fig_all.show()

    with open('C:\\Users\\RyanTrottier\\PycharmProjects\\untitled\\diffusion_barriers\\materials.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        for line in csv_dict.get_list():
            writer.writerow(line)
#%%
def make_transition(start, final, ts, t1, t2, output1, output2):
    from ase.io import write
    start = AseAtomsAdaptor.get_atoms(Poscar.from_dict(start['poscar']).structure)
    final = AseAtomsAdaptor.get_atoms(Poscar.from_dict(final['poscar']).structure)
    t1 = [AseAtomsAdaptor.get_atoms(Poscar.from_dict(x).structure) for x in t1['poscars']]
    t2 = [AseAtomsAdaptor.get_atoms(Poscar.from_dict(x).structure) for x in t2['poscars']]
    t1.reverse()
    ts = AseAtomsAdaptor.get_atoms(Poscar.from_dict(ts['poscar']).structure)
    write(output1, ts)
    write(output2, [start] + t1 + [ts] + t2 + [final])
    return

def write_transition(mat_dict, ts_dict, material, output_folder='D:\\Users\\RyanTrottier\\Documents\\Scrap\\paths'):
    energies = [x['energy'] for x in ts_dict[material]['lb']]
    i=energies.index(min(energies))
    mat = mat_dict[material]
    ts = ts_dict[material]
    make_transition(mat['start'], mat['finals'][i], ts['lb'][i], ts['ub'][i][0], ts['ub'][i][1],
                    os.path.join(output_folder, '{}.ts.cif'.format(material)), os.path.join(output_folder, '{}.path.xyz'.format(material)))
# for material in ['bafeo3', 'kbio3', 'mgal2o4', 'bahfo3', 'snhfo3', 'kbio3', 'zncro3']:
#     print(material)
#     write_transition(mat_dict, ts_dict, material)

def fermi_shift(material, i):
    base_vr = get_vasprun(mat_dict[material]['base']) # type: Vasprun
    start_vr = get_vasprun(mat_dict[material]['start']) # type: Vasprun
    ts_vr = get_vasprun(ts_dict[material]['lb'][i-1])
    return (base_vr.get_band_structure().get_band_gap()['energy'], ts_vr.efermi - base_vr.efermi, ts_vr.efermi - start_vr.efermi)

to_check = [
    ('laalo3', 2),
    ('znal2o4', 2),
    ('geni2o4', 1),
    ('mgal2o4', 1), # check
    ('scalo3', 1), # check
]