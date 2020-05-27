#%% Get Vaspruns
import numpy as np
import os
import Make_Dos
from Database_Tools import get_file

try:
    vaspruns
except:
    vaspruns = {}

materials_formated = [
'SnTaO3', 'GeNi2O4', 'MgAl2O4', 'InAlO3', 'CoAl2O4', 'ZnAl2O4', 'LaAlO3',
]

materials = [mat.lower() for mat in materials_formated]

i_list = [2, 1, 2, 4, 3, 3, 2,]

materials = zip(materials, i_list)

metals = ['bimno3']

for material, i in materials:
    print(material)
    # vaspruns[material] = {}
    for type in ['base', 'start', 'lb']:
        try:
            vaspruns[material][type]
        except:
            try:
                vaspruns[material]
            except:
                vaspruns[material] = {}
            if type == 'lb':
                v_file = get_file(fs, ts_dict[material][type][i-1]['vasprun'])
            else:
                v_file = get_file(fs, mat_dict[material][type]['vasprun'])
            vaspruns[material][type] = Vasprun(v_file)
            os.remove(v_file)


#%% Plot
import Make_Dos
from scipy.interpolate import spline
figsize=(3,5.5)
ylim = [0,3]
xlim = [-10,10]

    #vasprun shift name
plots = [



    # Fast
    # (vaspruns['bamn3po3']['base'], 0, 'BaMnO$_3$'),
    # (vaspruns['bamn3po3']['start'], 0, 'BaMnO$_3$ Vacancy'),
    # (vaspruns['bivo3']['base'], 0, 'BiVO3'),
    # (vaspruns['bivo3']['start'], 0, 'BiVO$_3$ Vacancy'),
    # (vaspruns['cual2o4']['base'], 0, 'CuAl2O4'),
    # (vaspruns['cual2o4']['start'], 0, 'CuAl$_2$O$_4$ Vacancy'),
    # (vaspruns['ktao3']['base'], 0, 'KTaO$_3$'),
    # (vaspruns['ktao3']['start'], 0, 'KTaO$_3$ Vacancy'),
    # (vaspruns['lasco3']['base'], 0, 'LaScO$_3$'),
    # (vaspruns['lasco3']['start'], 0, 'LaScO$_3$ Vacancy'),

    # Fast Metals
    # (vaspruns['batao3']['base'], 0, 'BaTaO$_3$'),
    # (vaspruns['batao3']['start'], 0, 'BaTaO$_3$ Vacancy'),
    # (vaspruns['bimno3']['base'], 0, 'BiMnO3'),
    # (vaspruns['bimno3']['start'], 0, 'BiMnO3 Vacancy'),
    # (vaspruns['lamn3po3']['base'], 0, 'LaMnO3'),
    # (vaspruns['lamn3po3']['start'], 0, 'LaMnO3 Vacancy'),
    # (vaspruns['mgcr2o4']['base'], -2, 'MgCr2O4'),
    # (vaspruns['mgcr2o4']['start'], 0, 'MgCr$_2$O$_4$ Vacancy'),
    # (vaspruns['camno3']['base'], 0, 'CaMnO$_3$'),
    # (vaspruns['camno3']['start'], 0, 'CaMnoO$_3$ Vacancy'),
    # (vaspruns['cacro3']['base'], 0, 'CaCrO$_3$'),
    # (vaspruns['cacro3']['start'], 0, 'CaCrO$_3$ Vacancy'),

    # Fast Needs Investigation
    # (vaspruns['baseo3']['base'], 0, 'BaSeO3'),
    # (vaspruns['baseo3']['start'], 0, 'BaSeO3 Vacancy'),
    # (vaspruns['nabio3']['base'], 0, 'NaBiO$_3$'),
    # (vaspruns['nabio3']['start'], 0, 'NaBiO$_3$ Vacancy'),




    # Slow
    # (vaspruns['bacoo3']['base'], 0, 'BaCoO$_3$'),
    # (vaspruns['bacoo3']['start'], 0, 'BaCoO$_3$ Vacancy'),
    # (vaspruns['bahfo3']['base'], 0, 'BaHfO$_3$'),
    # (vaspruns['bahfo3']['start'], 0, 'BaHfO$_3$ Vacancy'),
    # (vaspruns['lagao3']['base'], 0, 'LaGaO$_3$'),
    # (vaspruns['lagao3']['start'], 0, 'LaGaO$_3$ Vacancy'),
    # (vaspruns['mgal2o4']['base'], 0, 'MgAl2O4'),
    # (vaspruns['mgal2o4']['start'], 2.6, 'MgAl2O4 Vacancy'),


    # Slow Metals
    # (vaspruns['kbio3']['base'], 0, 'KBiO3'),
    # (vaspruns['kbio3']['start'], 0, 'KBiO3 Vacancy'),
    # (vaspruns['vmg2o4']['base'], 0, 'VMg$_2$O$_4$'),
    # (vaspruns['vmg2o4']['start'], 0, 'VMg$_2$O$_4$ Vacancy'),
    # (vaspruns['srcoo3']['base'], 0, 'SrCoO$_3$'),
    # (vaspruns['srcoo3']['start'], 0, 'SrCoO$_3$ Vacancy'),
    # (vaspruns['mnal2o4']['base'], 0, 'MnAl$_2$O$_4$'),
    # (vaspruns['mnal2o4']['start'], 0, 'MnAl$_2$O$_4$ Vacancy'),
    # (vaspruns['bafeo3']['base'], 0, 'BaFeO$_3$'),
    # (vaspruns['bafeo3']['start'], 0, 'BaFeO$_3$ Vacancy'),
    # (vaspruns['ynio3']['base'], 0, 'YNiO$_3$'),
    # (vaspruns['ynio3']['start'], 0, 'YNiO$_3$ Vacancy'),

    # ?
    # (vaspruns['basno3']['base'], 0, 'BaSnO3'),
    # (vaspruns['basno3']['start'], 0, 'BaSnO3 Vacancy'),

    # Retest
]

points = 1000
font = {'family' : 'Arial'}
colors = ['black', 'green', 'red']
# maximum = len(columns[1:])
maximum = 2

for vasprun, shift, name in plots:
    (title, columns, scaling_factors) = Make_Dos.make_dos(vasprun, ['O'])
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor('white')

    title_font = copy.deepcopy(font)
    title_font['size'] = 16
    for i in range(1,maximum,2):
        ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
        dos = np.array(columns[i]) - np.array(columns[i+1])
        energy = shift+np.array(columns[0])
        energy_new = np.linspace(energy.min(), energy.max(), points)
        # spl = make_interp_spline(x,y,k=3)
        dos_new = spline(energy,dos,energy_new)
        ax.plot(dos_new, energy_new, color=colors[int((i-1)/2)])
    ax.plot(ylim, [shift,shift], color='k', linestyle='--')
    # ax.legend([x[:-2] for x in title[1::2]] + ['Fermi'], loc='upper left')

    # Configure Axis
    label_font = copy.deepcopy(font)
    label_font['size'] = 14

    ax.set_ylim(xlim)
    ax.set_xlim(ylim)
    # ax.set_title(name, fontdict=title_font)
    ax.set_xticks([])
    # ax.set_yticks([])
    ax.yaxis.set_ticks_position('left')
    # ax.set_xlabel('DOS (arb. u.)', fontdict=label_font)
    # ax.set_ylabel('Energy (eV)', fontdict=label_font)

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

# for vasprun, shift, name in plots:
#     (title, columns, scaling_factors) = Make_Dos.make_dos(vasprun, [['O']])
#     fig, ax = plt.subplots(figsize=figsize)
#     for i in range(1,len(columns[1:]),2):
#         ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
#     ax.plot([shift,shift], ylim, color='k', linestyle='--')
#     ax.legend([x[:-2] for x in title[1::2]] + ['Fermi'], loc='upper left')
#     ax.set_ylim(ylim)
#     ax.set_xlim(xlim)
#     ax.set_title(name)
#     plt.show()


#%% DoS Overlay


plots = [
    # (vaspruns['geni2o4']['base'], vaspruns['geni2o4']['start'],  vaspruns['geni2o4']['lb'],
    #  0, 0, 'GeNi$_2$O$_4$'),
    (vaspruns['mgal2o4']['base'], vaspruns['mgal2o4']['start'],  vaspruns['mgal2o4']['lb'],
     0, 0, 'MgAl$_2$O$_4$'),
    (vaspruns['inalo3']['base'], vaspruns['inalo3']['start'],  vaspruns['inalo3']['lb'],
     0, 0, 'InAlO$_3$'),
    (vaspruns['coal2o4']['base'], vaspruns['coal2o4']['start'],  vaspruns['coal2o4']['lb'],
     0.2, 0.25, 'CoAl$_2$O$_4$'),
    (vaspruns['znal2o4']['base'], vaspruns['znal2o4']['start'],  vaspruns['znal2o4']['lb'],
     0, -0.2, 'ZnAl$_2$O$_4$'),
    (vaspruns['laalo3']['base'], vaspruns['laalo3']['start'],  vaspruns['laalo3']['lb'],
     0, 0, 'LaAlO$_3$'),
    # (vaspruns['sntao3']['base'], vaspruns['sntao3']['start'],  vaspruns['sntao3']['lb'],
    #  0, 0, 'SnTaO$_3$'),



]



import Make_Dos
from scipy.interpolate import spline
figsize=(4,6)
ylim = [0,1]
xlim = [-1,5]


    #vasprun shift name



points = 1000
font = {'family' : 'Arial'}
colors = ['black', 'green', 'red']
maximum = 2

for vasprun1, vasprun2, vasprun3, shift, shift2, name in plots:
    xlim = [0, vasprun1.get_band_structure().get_band_gap()['energy']]
    (title, columns, scaling_factors) = Make_Dos.make_dos(vasprun1, ['O'])
    (title2, columns2, scaling_factors2) = Make_Dos.make_dos(vasprun2, ['O'])
    (title3, columns3, scaling_factors3) = Make_Dos.make_dos(vasprun3, ['O'])
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor('white')

    title_font = copy.deepcopy(font)
    title_font['size'] = 16
    for i in range(1,maximum,2):
        # ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
        dos = np.array(columns[i]) - np.array(columns[i+1])
        energy = np.array(columns[0])
        energy_new = np.linspace(energy.min(), energy.max(), points)
        # spl = make_interp_spline(x,y,k=3)
        dos_new = spline(energy,dos,energy_new)
        ax.plot(dos_new, energy_new, color=colors[0])

        # ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
        dos2 = np.array(columns2[i]) - np.array(columns2[i+1])
        efermi2 = shift + (vasprun2.efermi - vasprun1.efermi)
        energy2 = shift+np.array(columns2[0])+(vasprun2.efermi - vasprun1.efermi)
        energy_new2 = np.linspace(energy2.min(), energy2.max(), points)
        # spl = make_interp_spline(x,y,k=3)
        dos_new2 = spline(energy2,dos2,energy_new2)
        ax.plot(dos_new2, energy_new2, color=colors[1], linestyle=':')


        # ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
        dos3 = np.array(columns3[i]) - np.array(columns3[i+1])
        efermi3 = shift2 + (vasprun3.efermi - vasprun1.efermi)
        energy3 = shift2+np.array(columns3[0])+(vasprun3.efermi - vasprun1.efermi)
        energy_new3 = np.linspace(energy3.min(), energy3.max(), points)
        # spl = make_interp_spline(x,y,k=3)
        dos_new3= spline(energy3,dos3,energy_new3)
        ax.plot(dos_new3, energy_new3, color=colors[2], linestyle='--')

    ax.legend(['Defect Free', 'O$_{Vac}$', 'TS'])
    # ax.plot(ylim, [0,0], color='k', linestyle='--')

    ax.plot(ylim, [0, 0], color=colors[0], linestyle='--')
    ax.plot(ylim, [efermi2, efermi2], color=colors[1], linestyle='--')
    ax.plot(ylim, [efermi3, efermi3], color=colors[2], linestyle='--')
    # ax.plot(ylim, [shift, shift], color='k', linestyle='--')
    # ax.legend([x[:-2] for x in title[1::2]] + ['Fermi'], loc='upper left')

    # Configure Axis
    label_font = copy.deepcopy(font)
    label_font['size'] = 14

    ax.set_ylim(xlim)
    ax.set_xlim(ylim)
    ax.set_title(name, fontdict=title_font)
    ax.set_xticks([])
    # ax.set_yticks([])
    ax.yaxis.set_ticks_position('left')
    # ax.set_xlabel('DOS (arb. u.)', fontdict=label_font)
    # ax.set_ylabel('Energy (eV)', fontdict=label_font)

    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)


    plt.show()
#%% Efermi
for material in materials:
    base = vaspruns[material]['base']
    vac = vaspruns[material]['start']
    print(material)
    print('Efermi {} -> {}'.format(base.efermi, vac.efermi))
    print('Band Gap {} -> {}'.format(base.get_band_structure().get_band_gap()['energy'], vac.get_band_structure().get_band_gap()['energy']))

# fig, ax = plt.subplots(figsize=(10,5))
#
# for column in columns[3:]:
#     ax.plot(columns[0], column)
# ax.legend(title[3:])
# ax.set_ylim([-1,1])
# ax.set_xlim([-10,10])
#
# maximum = np.argmax(np.array(columns[0])>0)
# # print(np.average(columns[0], weights=(np.array(columns[1])+np.array(columns[2]))[:maximum]) )
# fig.show()