

import scipy
cutoff = 5
min_cutoff = -100

x=np.array(csv_dict.get_column('ovac_final'))
y=np.array(csv_dict.get_column('lb'))/2+np.array(csv_dict.get_column('ub'))/2
x=x[np.where(y<cutoff)]
y=y[np.where(y<cutoff)]
linreg = scipy.stats.linregress(x, y)
prediction = linreg.intercept + linreg.slope * np.array(csv_dict.get_column('ovac_final'))


x=np.array(csv_dict.get_column('ts_efermi'))-np.array(csv_dict.get_column('efermi'))
# x= abs(x)
y=np.array(csv_dict.get_column('lb'))/2+np.array(csv_dict.get_column('ub'))/2
s=abs(np.array(csv_dict.get_column('ub'))-np.array(csv_dict.get_column('lb')))*500
marker='|'
plt.xlabel('E$_{fermi}$ Shift (eV)')
plt.ylabel('Diffusion Barrier (eV)')
plt.show()

# x=abs(np.array(csv_dict.get_column('formation_energy'))) + \
#   0.6*(np.array((csv_dict.get_column('o_p_band_center')) + 1.5/2*np.array(csv_dict.get_column('band_gap')))) + \
#   2.6*np.array(csv_dict.get_column('pauling_diff'))

# plt.ylabel('V$_{O}$ Energy (eV)')
# plt.xlabel('Deml')

# plt.ylabel('E$_{fermi}$ Shift (eV)')
# plt.xlabel('Band Gap (eV)')
#
# color=np.array(csv_dict.get_column('band_gap'))

# Band Gap
# np.array(csv_dict.get_column('band_gap'))

# Ovac
# np.array(csv_dict.get_column('ovac_final'))

# TS
# np.array(csv_dict.get_column('lb'))/2+np.array(csv_dict.get_column('ub'))/2
# prediction - ( np.array(csv_dict.get_column('lb'))/2+np.array(csv_dict.get_column('ub'))/2 )

# Fermi Shift
# np.array(csv_dict.get_column('ts_efermi'))-np.array(csv_dict.get_column('efermi'))
# np.array(csv_dict.get_column('efermi_vac'))-np.array(csv_dict.get_column('efermi'))
# np.array(csv_dict.get_column('ts_efermi'))-np.array(csv_dict.get_column('efermi_vac'))
#%%
x=np.array(csv_dict.get_column('ts_efermi'))-np.array(csv_dict.get_column('efermi_vac'))
y= prediction - ( np.array(csv_dict.get_column('lb'))/2+np.array(csv_dict.get_column('ub'))/2 )
color=np.array(csv_dict.get_column('band_gap'))
s=abs(np.array(csv_dict.get_column('ub'))-np.array(csv_dict.get_column('lb')))*500

cutoff = 5
min_cutoff = -5
s=s[np.where(y<cutoff)]
x=x[np.where(y<cutoff)]
color=color[np.where(y<cutoff)]
y=y[np.where(y<cutoff)]

s=s[np.where(y>min_cutoff)]
x=x[np.where(y>min_cutoff)]
color=color[np.where(y>min_cutoff)]
y=y[np.where(y>min_cutoff)]

color_cutoff = 0.5
s=s[np.where(color>color_cutoff)]
x=x[np.where(color>color_cutoff)]
y=y[np.where(color>color_cutoff)]
color=color[np.where(color>color_cutoff)]

s=5

z = np.polyfit(x, y, 1)
p = np.poly1d(z)
scat = plt.scatter(x, y, c=color, s=s, cmap='plasma', marker=marker)
plt.colorbar(scat, label='V$_O$ Energy (eV)')
plt.plot([min(x), max(x)], [p(min(x)), p(max(x))], "--", color='k')
plt.figure(figsize=(3,3))
plt.show()
print(scipy.stats.linregress(x, y))

# for property in csv_dict.get_header()[6:]:
#     y = prediction - (np.array(csv_dict.get_column('lb')) / 2 + np.array(csv_dict.get_column('ub')) / 2)
#     x = np.array(csv_dict.get_column(property))
#     x=x[np.where(y>min_cutoff)]
#     y=y[np.where(y>min_cutoff)]
#     linreg = scipy.stats.linregress(x, y)
#
#     print(property)
#     print(linreg.rvalue)
