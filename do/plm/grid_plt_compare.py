from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

#----------------------------------------------------------[Compare]

'''
str_tmp = 'dist_vp_o_1_2_fes_pka'
plm.grid_plt(
    list2d_file = [
        [f'../../../08.cn_o_h/330K/reweight/{str_tmp}.csv', 'Biasfactor 5.0'],
        [f'{str_tmp}.csv', 'Biasfactor 10.0'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = 'pKa',
    str_save = f'{str_tmp}.compare.pdf'
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes_deltag'
plm.grid_plt(
    list2d_file = [
        [f'../../../08.cn_o_h/330K/reweight/{str_tmp}.csv', 'Biasfactor 5.0'],
        [f'{str_tmp}.csv', 'Biasfactor 10.0'],
    ],
    str_xlabel = 'Time (ns)',
    str_ylabel = r'$\Delta$G (kJ/mol)',
    str_save = f'{str_tmp}.compare.pdf'
)
#'''
#'''
str_tmp = 'dhx_o_0_h_fes'
plm.grid_plt(
    list2d_file = [
        [f'../../../07.conformer/330K.2/reweight/{str_tmp}.grid', 'Biasfactor 5.0'],
        [f'{str_tmp}.grid', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'dhx(O$_0$H)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minzero = True,
    str_ylabel = 'FES (kJ/mol)',
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
plm.grid_plt(
    list2d_file = [
        [f'../../../08.cn_o_h/330K/reweight/{str_tmp}.grid', 'Biasfactor 5.0'],
        [f'{str_tmp}.grid', 'Biasfactor 10.0'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minzero = True,
    str_ylabel = 'FES (kJ/mol)',
    tup_xlim = (0,14),
    #tup_ylim = None,
)
#'''
'''
str_tmp = 'cn_o_h_fes'
plm.grid_plt(
    list2d_file = [
        [f'../../../10.cn_o_h/330K.0/sum_hills/{str_tmp}.5.dat', 'Biasfactor 5'],
        [f'{str_tmp}.5.dat', 'Biasfactor 10'],
    ],
    str_xlabel = r'cn(O$_C$H)',
    str_save = f'{str_tmp}.compare.pdf',
    bool_minus = True,
    bool_maxzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = None
)
#'''

plt.show()
