from pymatgen.io.vasp import Poscar, Vasprun
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen import Element
from get_migration import get_center_i, get_vacancy_diffusion_pathways_from_cell
import copy

# coefficients_001: -0.2088636566E+01
# 0.3434782476E+00 - 0.1303360998E+01
# 0.2474113532E+01
# Intercept_001: -0.5367985893E+01

def get_feat_1(features: dict):
    '''

    :param features:
    :return: float
    '''
# [(bond_A_length - unnormmed_diffusion_distance)]
    c = -0.2088636566e1
    bond_A_length = get_feat('bond_A_length', features)
    unnormmed_diffusion_distance = get_feat('diffusion_distance', features, False)
    return c*(bond_A_length - unnormmed_diffusion_distance)

def get_feat_2(features: dict):
    '''

    :param features:
    :return: float
    '''
# [(is_uncharged * unnormmed_band_gap)]
    c = 0.343478247
    is_uncharged = get_feat('is_uncharged', features, False)
    unnormmed_band_gap = get_feat('band_gap', features, False)
    return c*(is_uncharged * unnormmed_band_gap)

def get_feat_3(features: dict):
    '''

    :param features:
    :return: float
    '''
# [(bond_B_length * unnormmed_diffusion_distance)]
    c = - 0.1303360998E+01
    bond_B_length = get_feat('bond_B_length', features)
    unnormmed_diffusion_distance = get_feat('diffusion_distance', features, False)
    return c*(bond_B_length * unnormmed_diffusion_distance)

def get_feat_4(features: dict):
    '''

    :param features:
    :return: float
    '''
# [(formation_energy + atomization_energy)]
    c = 0.2474113532E+01
    formation_energy = get_feat('formation_energy', features)
    atomization_energy = get_feat('atomization_energy', features)
    return c*(formation_energy + atomization_energy)

chem_pots = {'H': -1.1217756,
 'Li': -1.90715197,
 'Be': -3.771544815,
 'C': -9.231242450625,
 'O': -4.94442311,
 'F': -1.2856885225,
 'Na': -1.338663779,
 'Mg': -1.5152837411111113,
 'Al': -3.74844693,
 'Si': -5.424965585,
 'P': -5.409410872857143,
 'S': -4.126883271875,
 'K': -1.0602799475,
 'Ca': -1.92335896,
 'Sc': -5.04495498,
 'Ti': -4.4520116666666665,
 'V': -4.26275489,
 'Cr': -5.94951935,
 'Mn': -6.921141819655173,
 'Fe': -5.18751316,
 'Co': -4.104854465,
 'Ni': -2.27729538,
 'Cu': -1.45330976,
 'Zn': -1.089421915,
 'Ga': -2.91921437,
 'Ge': -4.51841983,
 'Se': -3.505557640625,
 'Sr': -1.63689309,
 'Y': -5.31910215,
 'Zr': -6.58225052,
 'Nb': -7.03852874,
 'Cd': -0.736825565,
 'In': -2.5518483599999997,
 'Sn': -3.84632506,
 'Sb': -4.147950705,
 'Ba': -1.90871169,
 'La': -4.89907303,
 'Ce': -4.68467903,
 'Pr': -5.357426865,
 'Eu': -9.81043707,
 'Gd': -13.965768385,
 'Lu': -4.415364465,
 'Hf': -8.071571385,
 'Ta': -9.16697649,
 'Bi': -3.907564865}

atomization_energies = {'H': -0.00290213,
 'Li': -0.29803311,
 'Be': -0.03839057,
 'C': -1.3704101,
 'O': -1.90789148,
 'F': -0.62474621,
 'Na': -0.01109227,
 'Mg': -0.00049853,
 'Al': -0.14313117,
 'Si': -0.62906164,
 'P': -0.02135111,
 'S': -0.70537265,
 'K': -0.17837826,
 'Ca': -0.00693727,
 'Sc': -1.83549298,
 'Ti': -1.87254139,
 'V': -2.2389668,
 'Cr': -4.64837968,
 'Mn': -4.6922531,
 'Fe': -2.34548731,
 'Co': -0.802624,
 'Ni': 0.8892822,
 'Cu': 1.5261699,
 'Zn': -0.01114151,
 'Ga': -0.15173613,
 'Ge': -0.0192824,
 'Se': -0.89394831,
 'Sr': -0.02840861,
 'Y': -1.8921745,
 'Zr': -1.323768,
 'Cd': -0.01418287,
 'In': -0.25184654,
 'Sn': -0.02986457,
 'Sb': -0.87592437,
 'Ba': -0.03184493,
 'La': -0.61651049,
 'Ce': -0.7365616,
 'Pr': -2.55308713,
 'Nd': -4.44174128,
 'Sm': -6.80643846,
 'Eu': -8.12423019,
 'Gd': -9.92976694,
 'Tb': -8.47585603,
 'Dy': -7.24110832,
 'Ho': -6.69784917,
 'Er': -4.9557928,
 'Tm': -4.27717282,
 'Yb': 0.2062883,
 'Lu': -0.21764099,
 'Hf': -2.5648182,
 'Ta': -2.28858526,
 'Bi': -0.02466804}

def FormationEnergy(poscar, vasprun:Vasprun, elements):
    unit_s = poscar.structure
    elemental_energies = 0
    for atom in unit_s:
        element = atom.species_string
        elemental_energies += chem_pots[element]
    return (vasprun.final_energy - elemental_energies) / len(unit_s)


def AtomizationEnergy(poscar, vasprun:Vasprun, elements):
    unit_s = poscar.structure
    elemental_energies = 0
    for atom in unit_s:
        element = atom.species_string
        elemental_energies += atomization_energies[element]
    return -(vasprun.final_energy - elemental_energies) / len(unit_s)

def VoronoiInfo(poscar, vasprun, elements):
    base_s = poscar.structure
    start_i = get_center_i(base_s, Element('O'), )
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
    bonds = {
        'A_Bonds' : weight[elements[0]],
        'A_length' : total_length[elements[0]] / weight[elements[0]],
        'B_Bonds' : weight[elements[1]],
        'B_length' : total_length[elements[1]] / weight[elements[1]],
        'avg_length' : (total_length[elements[0]] + total_length[elements[1]]) / (weight[elements[0]] + weight[elements[1]]),
    }
    return bonds

def DiffusionDistance(poscar, vasprun, elements):

    base_s = poscar.structure

    start_i = get_center_i(base_s, Element('O'))
    start_coord = base_s.frac_coords[start_i]
    finals = get_vacancy_diffusion_pathways_from_cell(base_s, start_i, get_midpoints=True)[1]

    base_s.lattice.get_distance_and_image(start_coord, finals[0])
    distances = [2*base_s.lattice.get_distance_and_image(start_coord, final)[0] for final in finals]
    # print(distances)
    return distances

def Bandgap(poscar, vasprun, elements):
    bg = vasprun.get_band_structure().get_band_gap()['energy']
    return bg

def get_feat(feat, features, normalize=True):
    value = features[feat]
    if not normalize:
        return value
    normalization = {}
    normalization['bond_A_length'] = (1.8041439076909747, 3.453376021227605)
    normalization['bond_B_length'] = (1.6387592735086143, 2.9887241006183083)
    normalization['formation_energy'] = (-3.532782558000001, -1.0739940940000003)
    normalization['atomization_energy'] = (3.8513157839999996, 6.853106444)

    min_value = normalization[feat][0]
    max_value = normalization[feat][1]

    return (value-min_value)/(max_value-min_value)

def make_features(features, poscar, vasprun, elements):
    properties = [
        ('bond', VoronoiInfo(poscar, vasprun, elements)),
        ('diffusion_distance', DiffusionDistance(poscar, vasprun, elements)),
        ('band_gap', Bandgap(poscar, vasprun, elements)),
        ('formation_energy', FormationEnergy(poscar, vasprun, elements)),
        ('atomization_energy', AtomizationEnergy(poscar, vasprun, elements)),
    ]
    for name, property in properties:
        for feature_set in features.copy():
            value = property
            if type(value) == dict:
                for k in value:
                    feature_set['{}_{}'.format(name, k)] = value[k]
            if type(value) == list:  ### Will break for more than one list
                features = [copy.deepcopy(feature_set) for _ in range(len(value))]
                for i in range(len(value)):
                    features[i][name] = value[i]

            else:
                feature_set[name] = value
    return features


def get_diffusion_barrier(poscar : Poscar, vasprun : Vasprun, elements, unchargedP=1):
    '''

    :param poscar: pymatgen poscar file.  Ideally a unit structure
    :param vasprun:  vasprun corresponding to poscar
    :param elements: Elements of the A site and B site in the form [A_site, B_site]. i.e. for BaCoO3 it should be ['Ba', 'Co'].  O is not included
    :param unchargedP: whether diffusion barriers should be for charged or neutral diffusion
    :return: all possible diffusion barriers in the material
    '''
    features = [{'is_uncharged' : unchargedP}]
    feats = [get_feat_1, get_feat_2, get_feat_3, get_feat_4]
    barriers = []
    features = make_features(features, poscar, vasprun, elements)
    for feature_set in features:
        barrier=-0.5367985893E+01
        for feat in feats:
            barrier += feat(feature_set)
        barriers.append(barrier)

    return barriers
