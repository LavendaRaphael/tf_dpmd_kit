from tf_dpmd_kit import analysis

analysis.error(
    list_file = [
        'reweight.dist_vp_o_1_2_fes_pka.csv',
        '../03.metad/reweight.dist_vp_o_1_2_fes_pka.csv'
    ],
    str_save = f'reweight.dist_vp_o_1_2_fes_pka.ave.csv'
)
analysis.error(
    list_file = [
        'reweight.dist_vp_o_1_2_fes_deltag.csv',
        '../03.metad/reweight.dist_vp_o_1_2_fes_deltag.csv'
    ],
    str_save = f'reweight.dist_vp_o_1_2_fes_deltag.ave.csv'
)
