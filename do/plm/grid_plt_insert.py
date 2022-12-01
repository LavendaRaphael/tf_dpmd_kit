from tf_dpmd_kit import plm,analysis
import matplotlib.pyplot as plt

#----------------------------------------------------------[Compare]

fig, ax = plm.grid_plt(
    dict_data = {
        'fes': 'dist_vp_o_1_2_fes.grid',
    },
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    bool_minzero = True,
    tup_xlim = (None, 13),
    tup_ylim = (None, 60),
    bool_legend = False,
)

plm.insert_img(
    fig,ax,
    dict_img = {
        '023719.png': (0.1, 0.1, 0.2, 0.3),
        '023832.png': (0.2, 0.6, 0.3, 0.4),
        '025097.png': (0.5, 0.2, 0.3, 0.4),
    },
    dict_arrow = {},
    str_save = 'dist_vp_o_1_2_fes.img.pdf'
)

'''
fig, ax = plm.grid_plt(
    dict_data = {
        'fes': 'dhx_o_h_fes.grid',
    },
    str_xlabel = 'dhx(OH)',
    str_ylabel = 'FES (kJ/mol)',
    bool_minzero = True,
    tup_xlim = (-3, 3),
    tup_ylim = (None, 40),
    bool_legend = False,
)

plm.insert_img(
    fig,ax,
    dict_img = {
        '020000.png': (0.12, 0.5, 0.2, 0.3),
        '130000.png': (0.4, 0.5, 0.2, 0.3),
        '010000.png': (0.68, 0.6, 0.2, 0.3),
    },
    dict_arrow = {},
    str_save = 'dhx_o_h_fes.img.pdf'
)
'''

plt.show()


plt.show()
