import matplotlib.pyplot as plt
from tf_dpmd_kit import plm
from tf_dpmd_kit import plot

plot.set_rcparam()
cm = 1/2.54

dict_label = {
    'dist_o_0_h': r'R(O$_0$H$_0$) (Å)',
    'dist_o_2_h': r'R(O$_2$H$_1$) (Å)',
}

#'''
list_header = [
    'dist_o_0_h', 'dist_o_2_h'
]
for str_tmp in list_header:
    if str_tmp in dict_label:
        str_label = dict_label[str_tmp]
    else:
        str_label = str_tmp
    fig, ax = plm.colvar_hist_plt(
        str_header = str_tmp,
        str_label = str_label,
        tup_timerange = (0,28),
        tup_xlim = (0.9, 1.3),
    )
    plot.save(
        fig,
        tup_size = (8.6*cm, 4*cm),
        str_save = f'{str_tmp}.prob.svg',
    )
    plot.save(
        fig,
        tup_size = (8.6*cm, 4*cm),
        str_save = f'{str_tmp}.prob.pdf',
    )
#'''

plt.show()
