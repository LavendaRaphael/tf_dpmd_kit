from tf_dpmd_kit import plm
import math

'''
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
    dict_coef = {
        (0,2): -1,
        (2,14): 2
    },
    str_save = str_tmp+'_m4',
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
    dict_coef = {
        (0.7,1.4): 2,
        (1.4,2.5): -1
    },
    str_save = 'test',
    str_in = 'plm.in',
    float_volume = 15.6793091675**3
)
#'''
#'''
str_tmp = 'dhx_o_h_fes'
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
        5: f'{str_tmp}.grid',
    },
    dict_coef = {
        (-math.pi, -math.pi/4): -1,
        (-math.pi/4, math.pi/4): 1
    },
    str_save = str_tmp+'_CT',
    str_in = 'plm.in',
)

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
        5: f'{str_tmp}.grid',
    },
    dict_coef = {
        (-math.pi, -math.pi/4): -1,
        (math.pi/4, math.pi): 1
    },
    str_save = str_tmp+'_CC',
    str_in = 'plm.in',
)

#'''
