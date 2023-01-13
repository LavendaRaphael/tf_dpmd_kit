from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

str_tmp = 'dist_vp_o_1_2_fes'

plot.plt_compare(
    dict_data = {
        # '1ns':  f'analysis.0.{str_tmp}.grid',
        # '2ns':  f'analysis.1.{str_tmp}.grid',
        # '3ns':  f'analysis.2.{str_tmp}.grid',
        # '4ns':  f'analysis.3.{str_tmp}.grid',
          '5ns':  f'analysis.4.{str_tmp}.grid',
        # '6ns':  f'analysis.5.{str_tmp}.grid',
        # '7ns':  f'analysis.6.{str_tmp}.grid',
        # '8ns':  f'analysis.7.{str_tmp}.grid',
        # '9ns':  f'analysis.8.{str_tmp}.grid',
         '10ns':  f'analysis.9.{str_tmp}.grid',
        #'11ns': f'analysis.10.{str_tmp}.grid',
        #'12ns': f'analysis.11.{str_tmp}.grid',
        #'13ns': f'analysis.12.{str_tmp}.grid',
        #'14ns': f'analysis.13.{str_tmp}.grid',
         '15ns': f'analysis.14.{str_tmp}.grid',
        #'16ns': f'analysis.15.{str_tmp}.grid',
        #'17ns': f'analysis.16.{str_tmp}.grid',
        #'18ns': f'analysis.17.{str_tmp}.grid',
        #'19ns': f'analysis.18.{str_tmp}.grid',
         '20ns': f'analysis.19.{str_tmp}.grid',
        #'21ns': f'analysis.20.{str_tmp}.grid',
        #'22ns': f'analysis.21.{str_tmp}.grid',
        #'23ns': f'analysis.22.{str_tmp}.grid',
        #'24ns': f'analysis.23.{str_tmp}.grid',
         '25ns':             f'{str_tmp}.grid',
    },
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.converge.pdf',
    str_ylabel = 'FES (kJ/mol)',
    bool_minzero = True,
    #tup_xlim = (1.1, 8),
    tup_ylim = (None, 60)
)
plt.show()
