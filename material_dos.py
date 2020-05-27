#%% Get Vaspruns
import numpy as np
import os
import Make_Dos
from Database_Tools import get_file

try:
    vaspruns
except:
    vaspruns = {}

materials = [
    # 'laalo3', 'basno3', 'feal2o4', 'mgal2o4', 'mgcr2o4', 'cual2o4', 'bamn3po3', 'lamn3po3', 'bivo3', 'bimno3',
    # 'srmn3po3', 'srcoo3', 'baseo3', 'lasco3', 'ktao3', 'batao3', 'vmg2o4', 'kbio3',
    # 'mnal2o4', 'bafeo3', 'bacoo3',
    # 'camno3', 'ynio3', 'lagao3', 'bahfo3',
    'nabio3', 'cacro3'
]
metals = ['bimno3']

for material in materials:
    print(material)
    # vaspruns[material] = {}
    for type in ['base', 'start']:
        try:
            vaspruns[material][type]
        except:
            try:
                vaspruns[material]
            except:
                vaspruns[material] = {}
            if type == 'lb':
                v_file = get_file(fs, mat_dict[material][type]['vasprun'])
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
        # ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
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
    #Fast
    # (vaspruns['bamn3po3']['base'], vaspruns['bamn3po3']['start'], .5, 'BaMnO$_3$'),
    # (vaspruns['bivo3']['base'], vaspruns['bivo3']['start'], 0, 'BiVO$_3$'),
    # (vaspruns['cual2o4']['base'], vaspruns['cual2o4']['start'], 1.3, 'CuAl$_2$O$_4$'),
    # (vaspruns['ktao3']['base'], vaspruns['ktao3']['start'], .1, 'KTaO$_3$'),
    # (vaspruns['lasco3']['base'], vaspruns['lasco3']['start'], 2.9, 'LaScO$_3$'),
    #Slow
    # (vaspruns['bacoo3']['base'], vaspruns['bacoo3']['start'], .5, 'BaCoO$_3$'),
    # (vaspruns['bahfo3']['base'], vaspruns['bahfo3']['start'], 3, 'BaHfO$_3$'),
    # (vaspruns['lagao3']['base'], vaspruns['lagao3']['start'], 1.8, 'LaGaO$_3$'),
    # (vaspruns['mgal2o4']['base'], vaspruns['mgal2o4']['start'], 2.7, 'MgAl$_2$O$_4$'),
    #
    # (vaspruns['baseo3']['base'], vaspruns['baseo3']['start'], 1.7, 'BaSeO$_3$'),
    (vaspruns['nabio3']['base'], vaspruns['nabio3']['start'], 1, 'NaBiO$_3$'),


]



import Make_Dos
from scipy.interpolate import spline
figsize=(2,3)
ylim = [0,3]
xlim = [-5,5]


    #vasprun shift name



points = 1000
font = {'family' : 'Arial'}
colors = ['black', 'green', 'red']
maximum = 2

for vasprun1, vasprun2, shift, name in plots:
    (title, columns, scaling_factors) = Make_Dos.make_dos(vasprun1, ['O'])
    (title2, columns2, scaling_factors2) = Make_Dos.make_dos(vasprun2, ['O'])
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
        ax.plot(dos_new, energy_new, color=colors[int((i-1)/2)])

        # ax.plot(shift+np.array(columns[0]), np.array(columns[i]) - np.array(columns[i+1]))
        dos2 = np.array(columns2[i]) - np.array(columns2[i+1])
        energy2 = shift+np.array(columns2[0])
        energy_new2 = np.linspace(energy2.min(), energy2.max(), points)
        # spl = make_interp_spline(x,y,k=3)
        dos_new2 = spline(energy2,dos2,energy_new2)
        ax.plot(dos_new2, energy_new2, color=colors[int((i-1)/2)], linestyle=':')


    ax.legend(['Defect Free', 'O$_{Vac}$'])
    ax.plot(ylim, [0,0], color='k', linestyle='--')
    # ax.plot(ylim, [shift, shift], color='k', linestyle='--')
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