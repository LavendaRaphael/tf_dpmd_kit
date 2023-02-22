from tf_dpmd_kit import plm
import math

str_tmp = 'dist_vp_o_1_2_fes'
plm.get_pka_time(
    dict_file = {
         1: f'analysis.0.{str_tmp}.grid',
         2: f'analysis.1.{str_tmp}.grid',
         3: f'analysis.2.{str_tmp}.grid',
         4: f'analysis.3.{str_tmp}.grid',
         5: f'analysis.4.{str_tmp}.grid',
         6: f'analysis.5.{str_tmp}.grid',
         7: f'analysis.6.{str_tmp}.grid',
         8: f'analysis.7.{str_tmp}.grid',
         9: f'analysis.8.{str_tmp}.grid',
        10: f'analysis.9.{str_tmp}.grid',
        11: f'analysis.10.{str_tmp}.grid',
        12: f'analysis.11.{str_tmp}.grid',
        13: f'analysis.12.{str_tmp}.grid',
        14: f'analysis.13.{str_tmp}.grid',
        15: f'analysis.14.{str_tmp}.grid',
        16: f'analysis.15.{str_tmp}.grid',
        17: f'analysis.16.{str_tmp}.grid',
        18: f'analysis.17.{str_tmp}.grid',
        19: f'analysis.18.{str_tmp}.grid',
        20: f'analysis.19.{str_tmp}.grid',
        21: f'analysis.20.{str_tmp}.grid',
        22: f'analysis.21.{str_tmp}.grid',
        23: f'analysis.22.{str_tmp}.grid',
        24: f'analysis.23.{str_tmp}.grid',
        25:             f'{str_tmp}.grid',
    },
    dict_coef = {
        (0,2): -1,
        (2,14): 1
    },
    str_save = 'test'+str_tmp,
    str_in = 'plm.in',
    float_volume = 15.6793091675**3
)
