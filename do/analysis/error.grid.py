from tf_dpmd_kit import analysis

analysis.error(
    list_file = [
        'dist_vp_o_1_2_fes.grid',
        '../../../03.metad/320K/reweight/dist_vp_o_1_2_fes.grid'
    ],
    str_save = f'dist_vp_o_1_2_fes.ave.grid'
)
