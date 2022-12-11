from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

str_tmp = 'dist_vp_o_1_2_fes'
str_xlabel = r'R(V$_P$O$_C$) (Å)'

'''
plm.grid_plt(
    list2d_file = [
        [f'{str_tmp}.0.dat', '1ns'],
        [f'{str_tmp}.1.dat', '2ns'],
        [f'{str_tmp}.2.dat', '3ns'],
        [f'{str_tmp}.3.dat', '4ns'],
        [f'{str_tmp}.4.dat', '5ns'],
        #[f'{str_tmp}.5.dat', '6ns'],
        #[f'{str_tmp}.6.dat', '7ns'],
        #[f'{str_tmp}.7.dat', '8ns'],
        #[f'{str_tmp}.8.dat', '9ns'],
        #[f'{str_tmp}.9.dat', '10ns'],
    ],
    str_xlabel = str_xlabel,
    str_save = f'{str_tmp}.nosft.pdf',
    bool_minus = True,
    #bool_minzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = None
)
#'''
#'''
plm.grid_plt(
    list2d_file = [
        #[f'{str_tmp}.0.dat', '1ns'],
        [f'{str_tmp}.1.dat', '2ns'],
        #[f'{str_tmp}.2.dat', '3ns'],
        [f'{str_tmp}.3.dat', '4ns'],
        #[f'{str_tmp}.4.dat', '5ns'],
        [f'{str_tmp}.5.dat', '6ns'],
        #[f'{str_tmp}.6.dat', '7ns'],
        [f'{str_tmp}.7.dat', '8ns'],
        #[f'{str_tmp}.8.dat', '9ns'],
        [f'{str_tmp}.9.dat', '10ns'],
        #[f'{str_tmp}.10.dat', '11ns'],
        #[f'{str_tmp}.11.dat', '12ns'],
        #[f'{str_tmp}.12.dat', '13ns'],
        #[f'{str_tmp}.13.dat', '14ns'],
        #[f'{str_tmp}.14.dat', '15ns'],
    ],
    str_xlabel = str_xlabel,
    str_save = f'{str_tmp}.pdf',
    bool_minus = True,
    bool_maxzero = True,
    str_ylabel = 'Hills (kJ/mol)',
    tup_xlim = None,
    #tup_ylim = (None, 70)
)
#'''
plt.show()
