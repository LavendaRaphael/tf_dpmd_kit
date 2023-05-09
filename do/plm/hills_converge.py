from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

str_tmp = 'dist_vp_o_1_2_fes'

fig, ax = plot.plt_compare(
    dict_data = {
        # '1ns': f'{str_tmp}.0.dat',
        # '2ns': f'{str_tmp}.1.dat',
          '3ns': f'{str_tmp}.2.dat',
        # '4ns': f'{str_tmp}.3.dat',
        # '5ns': f'{str_tmp}.4.dat',
          '6ns': f'{str_tmp}.5.dat',
        # '7ns': f'{str_tmp}.6.dat',
        # '8ns': f'{str_tmp}.7.dat',
          '9ns': f'{str_tmp}.8.dat',
        #'10ns': f'{str_tmp}.9.dat',
        #'11ns': f'{str_tmp}.10.dat',
         '12ns': f'{str_tmp}.11.dat',
        #'13ns': f'{str_tmp}.12.dat',
        #'14ns': f'{str_tmp}.13.dat',
         '15ns': f'{str_tmp}.14.dat',
    },
    str_xlabel = r'R(V$_P$O$_{CA}$) (Å)',
    str_ylabel = 'V (kJ/mol)',
    bool_maxzero = True,
    bool_minus = True,
    #float_lw = 1,
)

plot.save(
    fig,
    tup_size = (8.6*cm, 5*cm),
    str_save = f'{str_tmp}.converge.svg',
)
plot.save(
    fig,
    tup_size = (8.6*cm, 5*cm),
    str_save = f'{str_tmp}.converge.pdf',
)

plt.show()
