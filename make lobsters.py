import os
os.environ['PMG_VASP_PSP_DIR'] = 'D:\\Users\\RyanTrottier\\Documents\\Scrap\\PMG'
os.environ['VASP_PSP_DIR'] = 'D:\\Users\\RyanTrottier\\Documents\\Scrap\\PMG'
import All_Materials
# from All_Materials import done,csv_dict,mat_dict
from ZachsMaterials import done, mat_dict
from AddDB import load_db
from Database_Tools import *
from Classes_Pymatgen import *
import pymongo
import numpy as np

#%%
l = 7
match_criteria = {
    # 'pathway_count': {'$exists' : False},
    # 'defect_type': {'$exists' : False},
    # 'ts_type': {'$exists' : False},
    'labels' : {'$all' : ['unit_cell'],
                '$nin' : ['surface']}
    # 'poscar.structure.lattice.a': {'$lt': l},
}

folder = 'D:\\Users\\RyanTrottier\\Documents\\Scrap\\lobsters'
(db,fs,client) = load_db()

# for material in csv_dict.materials.keys():
for material in mat_dict.keys():
    material = material.lower()
# for material in ['bicoo3']:
#     match_criteria['material'] = {'$all': [material], '$nin' : ['from_zach']}
    match_criteria['material'] = {'$all': [material] + ['from_zach']}
    runs = list(db.database.find(match_criteria).sort('energy', pymongo.ASCENDING))

    if 'ICOHPLIST_lobster' in runs[0]:
        continue
    print("{}: {}".format(material, len(runs)))
    if len(runs) > 0:
        [print(x['energy']) for x in runs]
    if len(runs) >= 1:
        run = runs[0]
        material_folder = os.path.join(folder, material)
        os.makedirs(material_folder, exist_ok=True)
        incar = Incar.from_dict(run['incar'])
        temp = get_file(fs, run['outcar'])
        magmom = [np.round(x['tot'],1) for x in Outcar(temp).magnetization]
        os.remove(temp)
        incar['MAGMOM'] = magmom
        incar['SYSTEM'] = material
        incar['KPAR'] = 3
        incar['NPAR'] = 3
        incar.write_file(os.path.join(material_folder, 'INCAR'))
        Kpoints.from_dict(run['kpoints']).write_file(os.path.join(material_folder, 'KPOINTS'))
        Poscar.from_dict(run['poscar']).write_file(os.path.join(material_folder, 'POSCAR'))
        Potcar(run['potcar']).write_file(os.path.join(material_folder, 'POTCAR'))
        with open(os.path.join(material_folder, 'DATABASE'), 'w') as f:
            f.write('''material {}
relaxation
unit_cell'''.format(' '.join(run['material'])))