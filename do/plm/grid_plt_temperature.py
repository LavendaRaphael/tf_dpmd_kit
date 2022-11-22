from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

#------------------------------------------------------------------[Temperature]

'''
#str_tmp = 'dist_vp_o_1_2_fes_deltag'
str_tmp = 'dhx_o_h_fes_CT_deltag'
plm.grid_plt(
    list2d_file = [
        [f'280K/reweight/{str_tmp}.csv', '280K'],
        [f'290K/reweight/{str_tmp}.csv', '290K'],
        [f'300K/reweight/{str_tmp}.csv', '300K'],
        [f'310K/reweight/{str_tmp}.csv', '310K'],
        [f'320K/reweight/{str_tmp}.csv', '320K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    str_save = f'reweight.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (-10, 70),
    tup_colormap = (280, 320)
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes_pka'
plm.grid_plt(
    list2d_file = [
        [f'280K/reweight/{str_tmp}.csv', '280K'],
        [f'290K/reweight/{str_tmp}.csv', '290K'],
        [f'300K/reweight/{str_tmp}.csv', '300K'],
        [f'310K/reweight/{str_tmp}.csv', '310K'],
        [f'320K/reweight/{str_tmp}.csv', '320K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'reweight.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 4.2),
    tup_colormap = (280, 320)
)
#'''
#'''
str_tmp='dist_vp_o_1_2_fes.10'
plm.grid_plt(
    list2d_file = [
        [f'280K.2/sum_hills/{str_tmp}.dat', '280K'],
        [f'290K.2/sum_hills/{str_tmp}.dat', '290K'],
        [f'300K.2/sum_hills/{str_tmp}.dat', '300K'],
        [f'310K.2/sum_hills/{str_tmp}.dat', '310K'],
        [f'320K.2/sum_hills/{str_tmp}.dat', '320K'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'Hills (kJ/mol)',
    str_save = f'sum_hills.{str_tmp}.pdf',
    bool_maxzero = True,
    bool_minus = True,
    tup_xlim = (0,14),
    tup_ylim = (None, None),
    tup_colormap = (280, 320)
)
#'''
'''
str_tmp = 'dhx_o_h_fes'
plm.grid_plt(
    list2d_file = [
        [f'280K/reweight/{str_tmp}.grid', '280K'],
        [f'290K/reweight/{str_tmp}.grid', '290K'],
        [f'300K/reweight/{str_tmp}.grid', '300K'],
        [f'310K/reweight/{str_tmp}.grid', '310K'],
        [f'320K/reweight/{str_tmp}.grid', '320K'],
    ],
    str_xlabel = r'dhx(O$_C$H)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight.{str_tmp}.pdf',
    bool_minzero = True,
    tup_colormap = (280, 320),
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
plm.grid_plt(
    list2d_file = [
        [f'280K/reweight/{str_tmp}.grid', '280K'],
        [f'290K/reweight/{str_tmp}.grid', '290K'],
        [f'300K/reweight/{str_tmp}.grid', '300K'],
        [f'310K/reweight/{str_tmp}.grid', '310K'],
        [f'320K/reweight/{str_tmp}.grid', '320K'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight.{str_tmp}.pdf',
    bool_minzero = True,
    tup_xlim = (0, 14),
    tup_ylim = (None, 60),
    tup_colormap = (280, 320),
)
#'''
plt.show()
