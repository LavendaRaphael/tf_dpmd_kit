from tf_dpmd_kit import plm

str_tmp = 'dist_vp_o_1_2_fes'
plm.get_pka_time(
    dict_file = {
        1: f'analysis.0.{str_tmp}.grid',
        2: f'analysis.1.{str_tmp}.grid',
        3: f'analysis.2.{str_tmp}.grid',
        4: f'analysis.3.{str_tmp}.grid',
        5: f'{str_tmp}.grid',
    },
    tup_xrange1 = (0,2),
    tup_xrange2 = (2,14),
    str_save = str_tmp,
    str_in = 'plm.in',
    float_volume = 15.6793091675**3
)

'''
run(
    str_file = 'dist_vp_c_fes.grid',
    tup_xrange1 = (1.16,1.66),
    tup_xrange2 = (2.56,'1M')
)

run(
    str_file = 'dist_vp_o_1_2_fes.grid',
    tup_xrange1 = (1.0,1.3),
    tup_xrange2 = (2.4,8.05)
)
'''
