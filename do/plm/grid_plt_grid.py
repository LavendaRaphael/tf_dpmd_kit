from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

#str_tmp = 'dist_vp_o_1_2_fes'
#str_xlabel = r'R(V$_P$O$_C$) (Å)'

#str_tmp = 'dhx_o_0_h_fes'
#str_xlabel = r'dhx(O$_0$H)'

#str_tmp = 'dhx_o_h_fes'
#str_xlabel = r'dhx(OH)'

str_tmp = 'cn_o_h_fes'
str_xlabel = r'cn(O$_C$H)'

plm.grid_plt(
    list2d_file = [
        [f'analysis.0.{str_tmp}.grid', '1ns'],
        [f'analysis.1.{str_tmp}.grid', '2ns'],
        [f'analysis.2.{str_tmp}.grid', '3ns'],
        [f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        #[f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        #[f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '5ns'],
    ],
    str_xlabel = str_xlabel,
    str_save = f'{str_tmp}.pdf',
    str_ylabel = 'FES (kJ/mol)',
    bool_minzero = True,
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
plt.show()
