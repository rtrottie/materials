#%% Imports

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
from pymatgen.analysis.local_env import VoronoiNN

from Make_Dos import make_dos
from material_helpers import *
try:
    client.close()
except:
    pass

(db,fs,client) = load_db()

#%% Load Materials

import All_Materials
import All_Charged_Materials
import ZachsMaterials
import ZungersMaterials

#%% Initialize CSV
csv_dict = Csv_dict()

add_no_ub = True

def addToCSV(mat_dict, ts_dict, unit_criteria, i_adj=0, extra_points={}):
    for material in list(mat_dict.keys())[:]:
        print(material)
        if material not in ['ceo2']:
            mat_unit_criteria = copy.deepcopy(unit_criteria)
            mat_unit_criteria['material']['$all'].append(material)
            print(mat_unit_criteria)
            units = list(db.database.find(mat_unit_criteria).sort('energy'))
            unit = units[0]
            mat = mat_dict[material]

            thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
            bulk = mat['base']
            start = mat['start']
            base = mat['base']
            finals = mat['finals']
            for i in range(1, len(mat['finals']) + 1):
                if thermo[i] and lb_bars[i] and ub_bars[i]:
                    ub = ts_dict[material]['ub'][i - 1]
                    ub_e = max(max([x['energies'] for x in ub])) - start['energy']
                elif thermo[i] and lb_bars[i]:
                    if not add_no_ub:
                        continue
                    ub = None
                    ub_e = None
                else:
                    continue
                lb = ts_dict[material]['lb'][i - 1]

                ovac_e = start['energy'] - bulk['energy']
                ovac_final = start['energy'] - bulk['energy']
                lb_e = lb['energy'] - start['energy']
                csv_dict.add_material(mat_dict[material]['base'], mat_dict[material]['start'],
                                      mat_dict[material]['finals'][i - 1], ovac_e, ovac_final, ub_e, lb_e, i+i_adj)
                for name, property in All_Materials.properties:
                    value = property.get(unit, base, start, finals[i - 1])
                    if type(value) == dict:
                        for k in value:
                            csv_dict.add_property(unit, '{}_{}'.format(name, k), value[k], i+i_adj)
                    else:
                        csv_dict.add_property(unit, name, value, i+i_adj)

                for name, property in All_Materials.ts_properties:
                    value = property.get(unit, base, start, finals[i - 1], lb, ub)
                    if type(value) == dict:
                        for k in value:
                            csv_dict.add_property(unit, '{}_{}'.format(name, k), value[k], i+i_adj)
                    else:
                        csv_dict.add_ts_property(unit, name, value, i+i_adj)

                for name, property in extra_points:
                    csv_dict.add_property(unit, name, value, i + i_adj)

#%% Add Main to CSV
addToCSV(All_Materials.mat_dict, All_Materials.ts_dict, All_Materials.unit_criteria)
print(len(csv_dict.get_list()))
#%% Make Charged to CSV
add_charged = True
if add_charged:
    addToCSV(All_Charged_Materials.mat_dict, All_Charged_Materials.ts_dict, All_Charged_Materials.unit_criteria, i_adj=100)
    print(len(csv_dict.get_list()))
#%% Make Zachs to CSV
addToCSV(ZachsMaterials.mat_dict, ZachsMaterials.ts_dict, ZachsMaterials.unit_criteria, i_adj=200)
print(len(csv_dict.get_list()))
#%% Make Zachs to CSV
addToCSV(ZungersMaterials.mat_dict, ZungersMaterials.ts_dict, ZungersMaterials.unit_criteria, i_adj=300)
print(len(csv_dict.get_list()))
#%% Load Charged Materials

#%% Write CSV

# with open('C:\\Users\\RyanTrottier\\PycharmProjects\\untitled\\diffusion_barriers\\materials.csv', 'w', newline='') as csvfile:
with open('C:\\Users\\RyanTrottier\\PycharmProjects\\materials\\data\\all_data\\materials_withubs.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    for line in csv_dict.get_list():
        writer.writerow(line)
