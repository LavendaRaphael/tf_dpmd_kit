from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

str_tmp = 'dist_vp_o_1_2_fes'
str_xlabel = r'R(V$_P$O$_C$) (Å)'

#str_tmp = 'dhx_o_h_fes'
#str_xlabel = r'dhx(OH)'

#str_tmp = 'cn_o_h_fes'
#str_xlabel = r'cn(O$_C$H)'

plm.grid_plt(
    dict_data = {
         '1ns':  f'analysis.0.{str_tmp}.grid',
         '2ns':  f'analysis.1.{str_tmp}.grid',
         '3ns':  f'analysis.2.{str_tmp}.grid',
         '4ns':  f'analysis.3.{str_tmp}.grid',
        # '5ns':  f'analysis.4.{str_tmp}.grid',
        # '6ns':  f'analysis.5.{str_tmp}.grid',
        # '7ns':  f'analysis.6.{str_tmp}.grid',
        # '8ns':  f'analysis.7.{str_tmp}.grid',
        # '9ns':  f'analysis.8.{str_tmp}.grid',
        #'10ns':  f'analysis.9.{str_tmp}.grid',
        #'11ns': f'analysis.10.{str_tmp}.grid',
        #'12ns': f'analysis.11.{str_tmp}.grid',
        #'13ns': f'analysis.12.{str_tmp}.grid',
        #'14ns': f'analysis.13.{str_tmp}.grid',
        '5ns':             f'{str_tmp}.grid',
    },
    str_xlabel = str_xlabel,
    str_save = f'{str_tmp}.pdf',
    str_ylabel = 'FES (kJ/mol)',
    bool_minzero = True,
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
plt.show()
