from tf_dpmd_kit import plm

#'''
str_tmp = 'dist_vp_o_1_2_fes'
plm.get_pka_time(
    dict_file = {
        1: f'analysis.0.{str_tmp}.grid',
        2: f'analysis.1.{str_tmp}.grid',
        3: f'analysis.2.{str_tmp}.grid',
        4: f'analysis.3.{str_tmp}.grid',
        #5: f'analysis.4.{str_tmp}.grid',
        #6: f'analysis.5.{str_tmp}.grid',
        #7: f'analysis.6.{str_tmp}.grid',
        #8: f'analysis.7.{str_tmp}.grid',
        #9: f'analysis.8.{str_tmp}.grid',
        #10: f'analysis.9.{str_tmp}.grid',
        #11: f'analysis.10.{str_tmp}.grid',
        #12: f'analysis.11.{str_tmp}.grid',
        #13: f'analysis.12.{str_tmp}.grid',
        #14: f'analysis.13.{str_tmp}.grid',
        5: f'{str_tmp}.grid',
    },
    tup_xrange1 = (0,2),
    tup_xrange2 = (2,14),
    str_save = str_tmp,
    str_in = 'plm.in',
    float_volume = 15.6793091675**3
)
#'''
'''
str_tmp = 'cn_o_h_fes'
plm.get_pka_time(
    dict_file = {
        1: f'analysis.0.{str_tmp}.grid',
        2: f'analysis.1.{str_tmp}.grid',
        3: f'analysis.2.{str_tmp}.grid',
        4: f'analysis.3.{str_tmp}.grid',
        5: f'{str_tmp}.grid',
    },
    tup_xrange1 = (1.4,2.5),
    tup_xrange2 = (0.7,1.4),
    str_save = str_tmp,
    str_in = 'plm.in',
    float_volume = 15.6793091675**3
)
'''
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
#'''
