from tf_dpmd_kit import plm

#------------------------------------------------------------------[Temperature]
'''
str_tmp = 'dist_vp_o_1_2_fes_deltag'
str_dir = 'sum_hills'
grid_plt(
    list2d_file = [
        [f'280K.1/{str_dir}/{str_tmp}.csv', '280K'],
        [f'290K.1/{str_dir}/{str_tmp}.csv', '290K'],
        [f'300K.1/{str_dir}/{str_tmp}.csv', '300K'],
        [f'310K.1/{str_dir}/{str_tmp}.csv', '310K'],
        [f'320K.1/{str_dir}/{str_tmp}.csv', '320K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    str_save = f'{str_dir}.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (-10, 70),
    tup_colormap = (280, 320)
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes_pka'
grid_plt(
    list2d_file = [
        [f'280K.1/{str_dir}/{str_tmp}.csv', '280K'],
        [f'290K.1/{str_dir}/{str_tmp}.csv', '290K'],
        [f'300K.1/{str_dir}/{str_tmp}.csv', '300K'],
        [f'310K.1/{str_dir}/{str_tmp}.csv', '310K'],
        [f'320K.1/{str_dir}/{str_tmp}.csv', '320K'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'{str_dir}.{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 4.2),
    tup_colormap = (280, 320)
)
#'''
'''
str_tmp='dist_vp_o_1_2_fes.10'
grid_plt(
    list2d_file = [
        [f'280K.1/sum_hills/{str_tmp}.dat', '280K'],
        [f'290K.1/sum_hills/{str_tmp}.dat', '290K'],
        [f'300K.1/sum_hills/{str_tmp}.dat', '300K'],
        [f'310K.1/sum_hills/{str_tmp}.dat', '310K'],
        [f'320K.1/sum_hills/{str_tmp}.dat', '320K'],
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
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        [f'290K/reweight_bias/{str_tmp}.grid', '290K'],
        [f'310K/reweight_bias/{str_tmp}.grid', '310K'],
        [f'330K/reweight_bias/{str_tmp}.grid', '330K'],
        [f'350K/reweight_bias/{str_tmp}.grid', '350K'],
        [f'370K/reweight_bias/{str_tmp}.grid', '370K'],
        [f'390K/reweight_bias/{str_tmp}.grid', '390K'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    bool_minzero = True,
    tup_xlim = (0, 14),
    tup_ylim = (None, 60),
    tup_colormap = (290, 390),
)
#'''
'''
str_tmp = 'dist_vp_c_fes'
#str_tmp = 'fes.dist_vp_c'
grid_plt(
    list2d_file = [
        [f'290K/reweight_bias/{str_tmp}.grid', '290K'],
        [f'310K/reweight_bias/{str_tmp}.grid', '310K'],
        [f'330K/reweight_bias/{str_tmp}.grid', '330K'],
        [f'350K/reweight_bias/{str_tmp}.grid', '350K'],
        [f'370K/reweight_bias/{str_tmp}.grid', '370K'],
        [f'390K/reweight_bias/{str_tmp}.grid', '390K'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_ylabel = 'FES (kJ/mol)',
    str_save = f'reweight_bias.{str_tmp}.pdf',
    tup_xlim = (0, 10),
    tup_ylim = (None, 60),
    bool_minzero = True,
    tup_colormap = (290, 390)
)
#'''
grid_plt(
    list2d_file = [
        [f'../../../03.390K_hills_reweight/390K/reweight_bias/{str_tmp}.grid', 'Biasfactor 5.0'],
        [f'{str_tmp}.grid', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minzero = True,
    str_ylabel = 'FES (kJ/mol)',
    tup_xlim = (0,10),
    tup_ylim = (None, 80),
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        [f'../../../03.390K_hills_reweight/390K/reweight_bias/{str_tmp}.grid', 'Biasfactor 5.0'],
        [f'{str_tmp}.grid', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minzero = True,
    str_ylabel = 'FES (kJ/mol)',
    tup_xlim = (0,14),
    tup_ylim = (None, 80),
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        ['../../../02.390K_hills/390K/sum_hills/fes.dist_vp_o_1_2.10.dat', 'Biasfactor 5.0'],
        [f'../../../04.hills_beta_10/390K/sum_hills/{str_tmp}.5.dat', 'Biasfactor 10.0'],
        [f'{str_tmp}.10.dat', 'MetaD'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minus = True,
    bool_maxzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = (0,14)
)
#'''
'''
str_tmp = 'dh_o_0_h_fes'
grid_plt(
    list2d_file = [
        [f'../../../07.conformer/330K/sum_hills/{str_tmp}.5.dat', 'Biasfactor 5'],
        [f'{str_tmp}.5.dat', 'Biasfactor 10'],
    ],
    str_xlabel = r'dh(O$_0$H)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minus = True,
    bool_maxzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = None
)
#'''
plt.show()
