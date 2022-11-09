from tf_dpmd_kit import plm

#----------------------------------------------------------[grid]
'''
str_tmp = 'dist_vp_c_hgx'
grid_plt(
    list2d_file = [
        #[f'analysis.0.{str_tmp}.grid', '1ns'],
        #[f'analysis.1.{str_tmp}.grid', '2ns'],
        #[f'analysis.2.{str_tmp}.grid', '3ns'],
        #[f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        [f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        [f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '10ns'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_ylabel = 'Probility Density',
    str_save = f'{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
'''
str_tmp = 'dist_vp_c_fes'
grid_plt(
    list2d_file = [
        #[f'analysis.0.{str_tmp}.grid', '1ns'],
        #[f'analysis.1.{str_tmp}.grid', '2ns'],
        #[f'analysis.2.{str_tmp}.grid', '3ns'],
        #[f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        [f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        [f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '10ns'],
    ],
    str_xlabel = r'R(V$_P$C) (Å)',
    str_save = f'{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
'''
str_tmp = 'dist_vp_o_1_2_fes'
grid_plt(
    list2d_file = [
        #[f'analysis.0.{str_tmp}.grid', '1ns'],
        #[f'analysis.1.{str_tmp}.grid', '2ns'],
        #[f'analysis.2.{str_tmp}.grid', '3ns'],
        #[f'analysis.3.{str_tmp}.grid', '4ns'],
        #[f'analysis.4.{str_tmp}.grid', '5ns'],
        [f'analysis.5.{str_tmp}.grid', '6ns'],
        #[f'analysis.6.{str_tmp}.grid', '7ns'],
        [f'analysis.7.{str_tmp}.grid', '8ns'],
        #[f'analysis.8.{str_tmp}.grid', '9ns'],
        [f'{str_tmp}.grid', '10ns'],
    ],
    str_xlabel = r'R(V$_P$O$_C$) (Å)',
    str_save = f'{str_tmp}.pdf',
    #tup_xlim = (1.1, 8),
    #tup_ylim = (None, 70)
)
#'''
plt.show()
