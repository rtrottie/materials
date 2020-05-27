from material_helpers import *
from Database_Tools import get_file
from All_Materials import *
import os
import shutil

folder = 'D:\\Users\\theme\\Documents\\scrap\\p2'
os.makedirs(folder, exist_ok=True)
o = 6

# for material in done:
for material in ['bialo3']:
    print(material)
    mat = mat_dict[material]

    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]

    ts_estimates = [gmean([lb_bars[i], ub_bars[i]]) if ub_bars[i] else 9999 for i in range(len(thermo))]
    lowest_ts_estimate = min(ts_estimates)
    i = ts_estimates.index(lowest_ts_estimate)-1
    i = 1
    if i == -1:
        continue

    base = mat['base']
    start = mat['start']
    final = mat['finals'][i]

    o_nelect = 6

    outcar_file = get_file(fs, base['outcar'])
    outcar = Outcar(outcar_file)
    os.remove(outcar_file)

    pairs = [
        ('base', base),
        ('start', start),
        ('final.{}'.format(i), final),
    ]

    material_folder = os.path.join(folder, material)
    for name, d in pairs:
        run_folder = os.path.join(material_folder, name)
        os.makedirs(run_folder, exist_ok=True)
        incar = Incar.from_dict(d['incar'])
        poscar = Poscar.from_dict(d['poscar'])
        kpoints = Kpoints.from_dict(d['kpoints'])
        outcar_file = get_file(fs, d['outcar'])
        outcar = Outcar(outcar_file)
        nelect = round(outcar.nelect-2)
        shutil.copy(outcar_file, os.path.join(run_folder, 'OUTCAR'))
        os.remove(outcar_file)

        incar['NELECT'] = nelect
        files = [
            ('INCAR', incar),
            ('KPOINTS', kpoints),
            ('POSCAR', poscar),
        ]
        for filename, pymg_object in files:
            pymg_object.write_file(os.path.join(run_folder, filename))