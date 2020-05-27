from material_helpers import *
from Database_Tools import get_file
from All_Charged_Materials import *
import os
import shutil

folder = 'D:\\Users\\theme\\Documents\\scrap\\p2_files'
os.makedirs(folder, exist_ok=True)
o = 6

for material in mat_dict:
    if mat_dict[material]['start'] == None:
        continue
    else:
        start = mat_dict[material]['start']
    final_exists = False
    for f in mat_dict[material]['finals']:
        if f != None:
            final_exists = True
            final = f
    if not final_exists:
        continue
    print(material)

    pairs = [
        ('start', start),
        ('final', final),
    ]

    material_folder = os.path.join(folder, material)
    for name, d in pairs:
        run_folder = os.path.join(material_folder, name)
        os.makedirs(run_folder, exist_ok=True)
        incar = Incar.from_dict(d['incar'])
        poscar = Poscar.from_dict(d['poscar'])
        kpoints = Kpoints.from_dict(d['kpoints'])
        for f in ['OUTCAR', 'vasprun']:
            outcar_file = get_file(fs, d[f.lower()])
            if f == 'vasprun':
                f = 'vasprun.xml'
            shutil.move(outcar_file, os.path.join(run_folder, f))

        files = [
            ('INCAR', incar),
            ('KPOINTS', kpoints),
            ('POSCAR', poscar),
        ]
        for filename, pymg_object in files:
            pymg_object.write_file(os.path.join(run_folder, filename))