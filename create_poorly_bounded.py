from material_helpers import *
from Database_Tools import get_file
from All_Materials import *
from Neb_Make import nebmake
import os
import shutil
import random
from Classes_Pymatgen import Incar,Poscar,Kpoints

folder = 'D:\\Users\\theme\\Documents\\scrap\\ts_benchmark'
os.makedirs(folder, exist_ok=True)
random.seed(1)
o = 6

cutoff = 2.46
max_range = 0.1

def write_run(run, run_folder, write_vasprun=False):
    os.makedirs(run_folder, exist_ok=True)
    incar = Incar.from_dict(run['incar'])
    poscar = Poscar.from_dict(run['poscar'])
    kpoints = Kpoints.from_dict(run['kpoints'])
    energy = run['energy']
    files = [
        ('INCAR', incar),
        ('KPOINTS', kpoints),
        ('POSCAR', poscar),
    ]
    for filename, pymg_object in files:
        pymg_object.write_file(os.path.join(run_folder, filename))
    with open(os.path.join(run_folder, 'energy.txt'), 'w') as f:
        f.write(str(energy))
    if write_vasprun:
        if not os.path.exists(os.path.join(run_folder, 'vasprun.xml')):
            print('  Writing Vasprun')
            get_file(fs, run['vasprun'], new_file=os.path.join(run_folder, 'vasprun.xml'))


for material in done:
# for material in []:
# for material in ['bialo3']:
    mat = mat_dict[material]
    ts = ts_dict[material]

    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]
    i = -1
    for lb_e, ub_e in zip(lb_bars, ub_bars):
        i += 1
        if i ==0:
            continue
        if lb_e and ub_e and random.random() < 0.2:
        # if (lb_e <= cutoff and ub_e >= cutoff) or (ub_e - lb_e > max_range and lb_e <= cutoff + max_range / 2 and ub_e >= cutoff - max_range / 2):
            print('{} {}'.format(material,i))
            material_folder = os.path.join(folder, material, str(i))
            lb = ts['lb'][i - 1]

            ubs = ts['ub'][i-1]
            start = mat['start']
            final = mat['finals'][i-1]

            #Write LB
            write_run(start, os.path.join(material_folder, 'start'), write_vasprun=True)
            write_run(final, os.path.join(material_folder, 'final'), write_vasprun=True)
            write_run(lb, os.path.join(material_folder, 'lb'), write_vasprun=True)
            write_run(mat['base'], os.path.join(material_folder, 'base'), write_vasprun=True)

            incar = Incar.from_dict(lb['incar'])
            incar['IMAGES'] = 1
            incar['SPRING'] = 0
            incar['LCLIMB'] = False
            incar['ISTART'] = 0
            incar['ICHARG'] = 2
            incar['IOPT'] = 1
            incar['IBRION'] = 3
            incar['POTIM'] = 0
            incar['EDIFFG'] = -0.05
            incar['EDIFF'] = 1e-5
            incar['NELM'] = 100
            incar['NSW'] = 5000

            kpoints = Kpoints.from_dict(lb['kpoints'])
            potcar = lb['potcar']
            relax = Poscar.from_dict(lb['poscar'])
            found_neb=False
            for tolerance in np.linspace(0.5, 2, 50):
                try:
                    neb = nebmake('', Poscar.from_dict(start['poscar']).structure,
                                  Poscar.from_dict(final['poscar']).structure, 2, tolerance=tolerance, write=False, quickfail=True)
                    found_neb=True
                    break
                except ValueError:
                    pass
            if not found_neb:
                print('Failed for {} {}'.format(material, i))
                continue

            start = Poscar(neb[0])
            relax = Poscar(neb[1])
            final = Poscar(neb[2])

            distances = [
                relax.structure.lattice.get_distance_and_image(relax.structure.frac_coords[i], start.structure.frac_coords[i])[0]
                for i in range(len(relax.structure))]
            diff_i = distances.index(max(distances))
            maxes = []
            distances = list(
                relax.structure.lattice.get_all_distances(relax.structure.frac_coords[diff_i], relax.structure.frac_coords)[0])
            maxes.append(distances.index(max(distances)))

            distances[maxes[-1]] = 0
            maxes.append(distances.index(max(distances)))

            sd = [[True] * 3 if i not in maxes else [False] * 3 for i in range(len(relax.structure))]
            relax = Poscar(relax.structure, selective_dynamics=sd)

            run_folder = os.path.join(material_folder, 'lb.neb')
            os.makedirs(run_folder, exist_ok=True)
            with open(os.path.join(run_folder, 'DATABASE'), 'w') as f:
                f.write('''material {}
relaxation
pathway_count {}
ts
ts_type pc neb
index {}'''.format(' '.join(lb['material']), lb['pathway_count'][0], lb['index'][0]))
            for f in ['', '00', '01', '02']:
                os.makedirs(os.path.join(run_folder,f), exist_ok=True)
            incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'lb',)
            files = [
                ('INCAR', incar),
                ('KPOINTS', kpoints),
                ('00/POSCAR', start),
                ('01/POSCAR', relax),
                ('02/POSCAR', final),
            ]
            fs_files = [
                (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
                (mat['start']['vasprun'],os.path.join(run_folder, '00', 'vasprun.xml')),
                (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '02', 'OUTCAR')),
                (mat['finals'][i-1]['vasprun'],os.path.join(run_folder, '02', 'vasprun.xml')),
            ]

            for oid, location in fs_files:
                if not os.path.exists(location):
                    get_file(fs, oid, new_file=location)

            for filename, pymg_object in files:
                pymg_object.write_file(os.path.join(run_folder, filename))
            # with open(os.path.join(run_folder, 'energy.txt'), 'w') as f:
            #     f.write(str(ub['energies'][ub_i]))
            with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
                f.write('\n'.join(potcar))

            neb = nebmake('', start.structure,
                                  final.structure, 4, tolerance=tolerance, write=False, quickfail=True)


            del incar['SPRING']
            incar['IMAGES'] = 3
            incar['LCLIMB'] = True

            run_folder = os.path.join(material_folder, 'neb.3')
            for f in ['', '00', '01', '02', '03', '04']:
                os.makedirs(os.path.join(run_folder,f), exist_ok=True)
            incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'neb')
            files = [
                ('INCAR', incar),
                ('KPOINTS', kpoints),
            ]


            fs_files = [
                (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
                (mat['start']['vasprun'],os.path.join(run_folder, '00', 'vasprun.xml')),
                (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '04', 'OUTCAR')),
                (mat['finals'][i-1]['vasprun'],os.path.join(run_folder, '04', 'vasprun.xml')),
            ]

            for oid, location in fs_files:
                if not os.path.exists(location):
                    get_file(fs, oid, new_file=location)


            for filename, pymg_object in files:
                pymg_object.write_file(os.path.join(run_folder, filename))

            with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
                f.write('\n'.join(potcar))
            for j, s in enumerate(neb):
                Poscar(s).write_file(os.path.join(run_folder, str(j).zfill(2), 'POSCAR'))





folder = 'D:\\Users\\theme\\Documents\\scrap\\fix_bounded'

for material in done:
# for material in []:
# for material in ['bialo3']:
    mat = mat_dict[material]
    ts = ts_dict[material]

    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]
    i = -1
    for lb_e, ub_e in zip(lb_bars, ub_bars):
        i += 1
        if i ==0:
            continue
        if (lb_e <= cutoff and ub_e >= cutoff):
        # if (lb_e <= cutoff and ub_e >= cutoff) or (ub_e - lb_e > max_range and lb_e <= cutoff + max_range / 2 and ub_e >= cutoff - max_range / 2):
            print('{} {}'.format(material,i))
            material_folder = os.path.join(folder, material, str(i))
            lb = ts['lb'][i - 1]
            ubs = ts['ub'][i-1]
            start = mat['start']
            final = mat['finals'][i-1]

            # Write LB
            write_run(start, os.path.join(material_folder, 'start'))
            write_run(final, os.path.join(material_folder, 'final'))
            write_run(lb, os.path.join(material_folder, 'lb'))
            j = 1
            ub_e1 = max(ubs[0]['energies'])
            ub_e2 = max(ubs[1]['energies'])
            if cutoff < ub_e - abs(ub_e1-ub_e2): # both upper bound above cutoff
                pass
            elif ub_e1 > ub_e2: # only ub1 is above cutoff
                ubs = ubs[:1]
            else:
                ubs = ubs[1:]
            for ub in ubs:
                ub_i = ub['energies'].index(max(ub['energies']))
                if ub_i == 0:
                    ub_i += 1
                elif ub_i == len(ub['energies'])-1:
                    ub_i -= 1
                incar = Incar.from_dict(lb['incar'])
                incar['IMAGES'] = 1
                incar['SPRING'] = 0
                kpoints = Kpoints.from_dict(lb['kpoints'])
                potcar = lb['potcar']
                init = Poscar.from_dict(ub['poscars'][ub_i-1])
                relax = Poscar.from_dict(ub['poscars'][ub_i])
                final = Poscar.from_dict(ub['poscars'][ub_i+1])

                distances = [
                    relax.structure.lattice.get_distance_and_image(relax.structure.frac_coords[i], init.structure.frac_coords[i])[0]
                    for i in range(len(relax.structure))]
                diff_i = distances.index(max(distances))
                maxes = []
                distances = list(
                    relax.structure.lattice.get_all_distances(relax.structure.frac_coords[diff_i], relax.structure.frac_coords)[0])
                maxes.append(distances.index(max(distances)))

                distances[maxes[-1]] = 0
                maxes.append(distances.index(max(distances)))

                sd = [[True] * 3 if i not in maxes else [False] * 3 for i in range(len(relax.structure))]
                relax = Poscar(relax.structure, selective_dynamics=sd)

                run_folder = os.path.join(material_folder, 'neb.{}'.format(ub['interpolation_direction'][0]))
                for f in ['', '00', '01', '02']:
                    os.makedirs(os.path.join(run_folder,f), exist_ok=True)
                incar['SYSTEM'] = '{} {} {} {}'.format(ub['material'][1], i-1, 'neb', ub['interpolation_direction'][0])
                files = [
                    ('INCAR', incar),
                    ('KPOINTS', kpoints),
                    ('00/POSCAR', init),
                    ('01/POSCAR', relax),
                    ('02/POSCAR', final),
                ]
                for filename, pymg_object in files:
                    pymg_object.write_file(os.path.join(run_folder, filename))
                with open(os.path.join(run_folder, 'energy.txt'), 'w') as f:
                    f.write(str(ub['energies'][ub_i]))
                with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
                    f.write('\n'.join(potcar))


            # Write Pathway


            # write example pathway

            # create recommended minimization
raise Exception('Done')
folder = 'D:\\Users\\theme\\Documents\\scrap\\redo_lb'
for material in done:
# for material in []:
# for material in ['bialo3']:
    mat = mat_dict[material]
    ts = ts_dict[material]

    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]
    i = -1
    for lb_e, ub_e in zip(lb_bars, ub_bars):
        i += 1
        if i ==0:
            continue
        if lb_e and ub_e and ub_e - lb_e > 0.25 and 'IMAGES' not in ts['lb'][i-1]['incar']:
        # if (lb_e <= cutoff and ub_e >= cutoff) or (ub_e - lb_e > max_range and lb_e <= cutoff + max_range / 2 and ub_e >= cutoff - max_range / 2):
            print('{} {}'.format(material,i))
            material_folder = os.path.join(folder, material, str(i))
            lb = ts['lb'][i - 1]

            ubs = ts['ub'][i-1]
            start = mat['start']
            final = mat['finals'][i-1]

            # Write LB
            # write_run(start, os.path.join(material_folder, 'start'), write_vasprun=True)
            # write_run(final, os.path.join(material_folder, 'final'), write_vasprun=True)
            # write_run(lb, os.path.join(material_folder, 'lb'), write_vasprun=True)
            # write_run(mat['base'], os.path.join(material_folder, 'base'), write_vasprun=True)

            # incar = Incar.from_dict(lb['incar'])
            # incar['IMAGES'] = 1
            # incar['SPRING'] = 0
            # kpoints = Kpoints.from_dict(lb['kpoints'])
            # potcar = lb['potcar']
            # relax = Poscar.from_dict(lb['poscar'])
            # found_neb=False
            # for tolerance in np.linspace(0.5, 2, 50):
            #     try:
            #         neb = nebmake('', Poscar.from_dict(start['poscar']).structure,
            #                       Poscar.from_dict(final['poscar']).structure, 2, tolerance=tolerance, write=False, quickfail=True)
            #         found_neb=True
            #         break
            #     except ValueError:
            #         pass
            # if not found_neb:
            #     print('Failed for {} {}'.format(material, i))
            #     continue
            #
            # start = Poscar(neb[0])
            # relax = Poscar(neb[1])
            # final = Poscar(neb[2])
            #
            # distances = [
            #     relax.structure.lattice.get_distance_and_image(relax.structure.frac_coords[i], start.structure.frac_coords[i])[0]
            #     for i in range(len(relax.structure))]
            # diff_i = distances.index(max(distances))
            # maxes = []
            # distances = list(
            #     relax.structure.lattice.get_all_distances(relax.structure.frac_coords[diff_i], relax.structure.frac_coords)[0])
            # maxes.append(distances.index(max(distances)))
            #
            # distances[maxes[-1]] = 0
            # maxes.append(distances.index(max(distances)))
            #
            # sd = [[True] * 3 if i not in maxes else [False] * 3 for i in range(len(relax.structure))]
            # relax = Poscar(relax.structure, selective_dynamics=sd)

            run_folder = os.path.join(material_folder, 'lb.neb')
            os.makedirs(run_folder, exist_ok=True)
            with open(os.path.join(run_folder, 'DATABASE'), 'w') as f:
                f.write('''material {}
relaxation
pathway_count {}
ts
ts_type pc neb
index {}'''.format(' '.join(lb['material']), lb['pathway_count'][0], lb['index'][0]))
            # for f in ['', '00', '01', '02']:
            #     os.makedirs(os.path.join(run_folder,f), exist_ok=True)
            # incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'lb',)
            # files = [
            #     ('INCAR', incar),
            #     ('KPOINTS', kpoints),
            #     ('00/POSCAR', start),
            #     ('01/POSCAR', relax),
            #     ('02/POSCAR', final),
            # ]
            # fs_files = [
            #     (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
            #     (mat['start']['vasprun'],os.path.join(run_folder, '00', 'vasprun.xml')),
            #     (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '02', 'OUTCAR')),
            #     (mat['finals'][i-1]['vasprun'],os.path.join(run_folder, '02', 'vasprun.xml')),
            # ]
            #
            # for oid, location in fs_files:
            #     if not os.path.exists(location):
            #         get_file(fs, oid, new_file=location)
            #
            # for filename, pymg_object in files:
            #     pymg_object.write_file(os.path.join(run_folder, filename))
            # # with open(os.path.join(run_folder, 'energy.txt'), 'w') as f:
            # #     f.write(str(ub['energies'][ub_i]))
            # with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
            #     f.write('\n'.join(potcar))
            #
            # neb = nebmake('', start.structure,
            #                       final.structure, 4, tolerance=tolerance, write=False, quickfail=True)
            #
            #
            # del incar['SPRING']
            # incar['IMAGES'] = 3
            # incar['LCLIMB'] = True
            #
            # run_folder = os.path.join(material_folder, 'neb.3')
            # for f in ['', '00', '01', '02', '03', '04']:
            #     os.makedirs(os.path.join(run_folder,f), exist_ok=True)
            # incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'neb')
            # files = [
            #     ('INCAR', incar),
            #     ('KPOINTS', kpoints),
            # ]
            #
            #
            # fs_files = [
            #     (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
            #     (mat['start']['vasprun'],os.path.join(run_folder, '00', 'vasprun.xml')),
            #     (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '04', 'OUTCAR')),
            #     (mat['finals'][i-1]['vasprun'],os.path.join(run_folder, '04', 'vasprun.xml')),
            # ]
            #
            # for oid, location in fs_files:
            #     if not os.path.exists(location):
            #         get_file(fs, oid, new_file=location)
            #
            #
            # for filename, pymg_object in files:
            #     pymg_object.write_file(os.path.join(run_folder, filename))
            #
            # with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
            #     f.write('\n'.join(potcar))
            # for j, s in enumerate(neb):
            #     Poscar(s).write_file(os.path.join(run_folder, str(j).zfill(2), 'POSCAR'))

folder = 'D:\\Users\\theme\\Documents\\scrap\\for_table'
for material in done:
    mat = mat_dict[material]
    ts = ts_dict[material]

    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]
    i = -1
    for lb_e, ub_e in zip(lb_bars, ub_bars):
        i += 1
        if i ==0:
            continue
        if lb_e and ub_e and ub_e - lb_e < 0.25 and 'IMAGES' not in ts['lb'][i-1]['incar']:
        # if (lb_e <= cutoff and ub_e >= cutoff) or (ub_e - lb_e > max_range and lb_e <= cutoff + max_range / 2 and ub_e >= cutoff - max_range / 2):
            print('{} {}'.format(material,i))
            material_folder = os.path.join(folder, material, str(i))
            lb = ts['lb'][i - 1]
            ubs = ts['ub'][i-1]
            start = mat['start']
            final = mat['finals'][i-1]

            # Write LB
            write_run(start, os.path.join(material_folder, 'start'), write_vasprun=False)
            write_run(final, os.path.join(material_folder, 'final'), write_vasprun=False)
            write_run(lb, os.path.join(material_folder, 'lb'), write_vasprun=False)
            write_run(mat['base'], os.path.join(material_folder, 'base'), write_vasprun=False)

            incar = Incar.from_dict(lb['incar'])
            incar['IMAGES'] = 1
            incar['SPRING'] = 0
            kpoints = Kpoints.from_dict(lb['kpoints'])
            potcar = lb['potcar']
            relax = Poscar.from_dict(lb['poscar'])
            found_neb=False
            for tolerance in np.linspace(0.5, 2, 50):
                try:
                    neb = nebmake('', Poscar.from_dict(start['poscar']).structure,
                                  Poscar.from_dict(final['poscar']).structure, 2, tolerance=tolerance, write=False, quickfail=True)
                    found_neb=True
                    break
                except ValueError:
                    pass
            if not found_neb:
                print('Failed for {} {}'.format(material, i))
                continue

            start = Poscar(neb[0])
            relax = Poscar(neb[1])
            final = Poscar(neb[2])

            distances = [
                relax.structure.lattice.get_distance_and_image(relax.structure.frac_coords[i], start.structure.frac_coords[i])[0]
                for i in range(len(relax.structure))]
            diff_i = distances.index(max(distances))
            maxes = []
            distances = list(
                relax.structure.lattice.get_all_distances(relax.structure.frac_coords[diff_i], relax.structure.frac_coords)[0])
            maxes.append(distances.index(max(distances)))

            distances[maxes[-1]] = 0
            maxes.append(distances.index(max(distances)))

            sd = [[True] * 3 if i not in maxes else [False] * 3 for i in range(len(relax.structure))]
            relax = Poscar(relax.structure, selective_dynamics=sd)

            run_folder = os.path.join(material_folder, 'lb.neb')
            for f in ['', '00', '01', '02']:
                os.makedirs(os.path.join(run_folder,f), exist_ok=True)
            incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'lb')
            files = [
                ('INCAR', incar),
                ('KPOINTS', kpoints),
                ('00/POSCAR', start),
                ('01/POSCAR', relax),
                ('02/POSCAR', final),
            ]

            fs_files = [
                (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
                (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '02', 'OUTCAR')),
            ]
            for oid, location in fs_files:
                if not os.path.exists(location):
                    get_file(fs, oid, new_file=location)
            for filename, pymg_object in files:
                pymg_object.write_file(os.path.join(run_folder, filename))
            with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
                f.write('\n'.join(potcar))

            neb = nebmake('', start.structure,
                                  final.structure, 4, tolerance=tolerance, write=False, quickfail=True)
            del incar['SPRING']
            incar['IMAGES'] = 3
            incar['LCLIMB'] = True
            incar['ISTART'] = 0
            incar['ICHARG'] = 2
            incar['NSW'] = 5000
            incar['IBRION'] = 3
            incar['IOPT'] = 1
            incar['POTIM'] = 0
            incar['EDIFFG'] = -0.05

            run_folder = os.path.join(material_folder, 'neb.3')
            for f in ['', '00', '01', '02', '03', '04']:
                os.makedirs(os.path.join(run_folder,f), exist_ok=True)
            incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'neb')
            files = [
                ('INCAR', incar),
                ('KPOINTS', kpoints),
            ]
            fs_files = [
                (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
                (mat['start']['vasprun'],os.path.join(run_folder, '00', 'vasprun.xml')),
                (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '04', 'OUTCAR')),
                (mat['finals'][i-1]['vasprun'],os.path.join(run_folder, '04', 'vasprun.xml')),
            ]

            for oid, location in fs_files:
                if not os.path.exists(location):
                    get_file(fs, oid, new_file=location)
            for filename, pymg_object in files:
                pymg_object.write_file(os.path.join(run_folder, filename))
            with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
                f.write('\n'.join(potcar))
            for j, s in enumerate(neb):
                Poscar(s).write_file(os.path.join(run_folder, str(j).zfill(2), 'POSCAR'))
#%%
folder = 'D:\\Users\\theme\\Documents\\scrap\\ts_dos'
predictions = [('mgal2o4', 2), ('inalo3', 4), ('coal2o4', 3), ('znal2o4', 3), ('laalo3', 2), ('catio3', 1), ('srtio3', 2), ('bialo3', 1), ('snhfo3', 3), ('mgseo3', 3)]

for material, i in predictions:
    mat = mat_dict[material]
    ts = ts_dict[material]

    thermo, lb_bars, ub_bars, delta = get_bars(mat, ts_dict[material])
    bars_list = [thermo, lb_bars, ub_bars]
        # if (lb_e <= cutoff and ub_e >= cutoff) or (ub_e - lb_e > max_range and lb_e <= cutoff + max_range / 2 and ub_e >= cutoff - max_range / 2):
    print('{} {}'.format(material,i))
    material_folder = os.path.join(folder, material, str(i))
    lb = ts['lb'][i - 1]
    ubs = ts['ub'][i-1]
    start = mat['start']
    final = mat['finals'][i-1]

    # Write LB
    write_run(start, os.path.join(material_folder, 'start'), write_vasprun=False)
    write_run(final, os.path.join(material_folder, 'final'), write_vasprun=False)
    write_run(lb, os.path.join(material_folder, 'lb'), write_vasprun=False)
    write_run(mat['base'], os.path.join(material_folder, 'base'), write_vasprun=False)

    incar = Incar.from_dict(lb['incar'])
    incar['IMAGES'] = 1
    incar['SPRING'] = 0
    kpoints = Kpoints.from_dict(lb['kpoints'])
    potcar = lb['potcar']
    relax = Poscar.from_dict(lb['poscar'])
    found_neb=False
    for tolerance in np.linspace(0.5, 2, 50):
        try:
            neb = nebmake('', Poscar.from_dict(start['poscar']).structure,
                          Poscar.from_dict(final['poscar']).structure, 2, tolerance=tolerance, write=False, quickfail=True)
            found_neb=True
            break
        except ValueError:
            pass
    if not found_neb:
        print('Failed for {} {}'.format(material, i))
        continue

    start = Poscar(neb[0])
    relax = Poscar(neb[1])
    final = Poscar(neb[2])

    distances = [
        relax.structure.lattice.get_distance_and_image(relax.structure.frac_coords[i], start.structure.frac_coords[i])[0]
        for i in range(len(relax.structure))]
    diff_i = distances.index(max(distances))
    maxes = []
    distances = list(
        relax.structure.lattice.get_all_distances(relax.structure.frac_coords[diff_i], relax.structure.frac_coords)[0])
    maxes.append(distances.index(max(distances)))

    distances[maxes[-1]] = 0
    maxes.append(distances.index(max(distances)))

    sd = [[True] * 3 if i not in maxes else [False] * 3 for i in range(len(relax.structure))]
    relax = Poscar(relax.structure, selective_dynamics=sd)

    run_folder = os.path.join(material_folder, 'lb.neb')
    for f in ['', '00', '01', '02']:
        os.makedirs(os.path.join(run_folder,f), exist_ok=True)
    incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'lb')
    files = [
        ('INCAR', incar),
        ('KPOINTS', kpoints),
        ('00/POSCAR', start),
        ('01/POSCAR', relax),
        ('02/POSCAR', final),
    ]

    fs_files = [
        (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
        (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '02', 'OUTCAR')),
    ]
    for oid, location in fs_files:
        if not os.path.exists(location):
            get_file(fs, oid, new_file=location)
    for filename, pymg_object in files:
        pymg_object.write_file(os.path.join(run_folder, filename))
    with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
        f.write('\n'.join(potcar))

    neb = nebmake('', start.structure,
                          final.structure, 4, tolerance=tolerance, write=False, quickfail=True)
    del incar['SPRING']
    incar['IMAGES'] = 3
    incar['LCLIMB'] = True
    incar['ISTART'] = 0
    incar['ICHARG'] = 2
    incar['NSW'] = 5000
    incar['IBRION'] = 3
    incar['IOPT'] = 1
    incar['POTIM'] = 0
    incar['EDIFFG'] = -0.05

    run_folder = os.path.join(material_folder, 'neb.3')
    for f in ['', '00', '01', '02', '03', '04']:
        os.makedirs(os.path.join(run_folder,f), exist_ok=True)
    incar['SYSTEM'] = '{} {} {}'.format(material, i-1, 'neb')
    files = [
        ('INCAR', incar),
        ('KPOINTS', kpoints),
    ]
    fs_files = [
        (mat['start']['outcar'], os.path.join(run_folder, '00', 'OUTCAR')),
        (mat['start']['vasprun'],os.path.join(run_folder, '00', 'vasprun.xml')),
        (mat['finals'][i-1]['outcar'],os.path.join(run_folder, '04', 'OUTCAR')),
        (mat['finals'][i-1]['vasprun'],os.path.join(run_folder, '04', 'vasprun.xml')),
    ]

    for oid, location in fs_files:
        if not os.path.exists(location):
            get_file(fs, oid, new_file=location)
    for filename, pymg_object in files:
        pymg_object.write_file(os.path.join(run_folder, filename))
    with open(os.path.join(run_folder, 'potcar.txt'), 'w') as f:
        f.write('\n'.join(potcar))
    for j, s in enumerate(neb):
        Poscar(s).write_file(os.path.join(run_folder, str(j).zfill(2), 'POSCAR'))

    with open(os.path.join(run_folder, 'DATABASE'), 'w') as f:
            f.write('''material {}
relaxation
pathway_count {}
ts
ts_type neb
index {}'''.format(' '.join(lb['material']), lb['pathway_count'][0], lb['index'][0]))