import matplotlib.pyplot as plt
from tf_dpmd_kit import plm
from tf_dpmd_kit import plot

plot.set_rcparam()
cm = 1/2.54

dict_label = {
    'dist_o_0_h': r'R(O$_0$H$_0$) (Å)',
    'dist_o_2_h': r'R(O$_2$H$_1$) (Å)',
    'dist_vp_o_1_2': r'R(V$_p$O$_{CA}$) (Å)',
    'dhx_o_0_h': r"α'(O$_2$-C-O$_0$-H$_0$) (rad)",
    'cnx_o_0_h': r"cn(O$_0$H$_W$)"
}

#'''
list_header = [
    'dist_o_0_h', 
    'dist_o_2_h',
    'dist_vp_o_1_2',
    'dhx_o_0_h',
    'cnx_o_0_h'
]
tup_size = (4.3*cm, 4*cm)

dict_xlim = {
    'dist_o_0_h': (0.9, 1.3),
    'dist_o_2_h': (0.9, 1.3),
    'dist_vp_o_1_2': (1, 1.4),
    'dhx_o_0_h': (None, 3),
    'cnx_o_0_h': (None, 0.03)
}

for str_tmp in list_header:
    if str_tmp in dict_label:
        str_label = dict_label[str_tmp]
    else:
        str_label = str_tmp
    if str_tmp in dict_xlim:
        tup_xlim = dict_xlim[str_tmp]
    else:
        tup_xlim = None
    fig, ax = plm.colvar_hist_plt(
        str_header = str_tmp,
        str_label = str_label,
        tup_timerange = (0,28),
        tup_xlim = tup_xlim,
    )
    plot.save(
        fig,
        tup_size = tup_size,
        str_save = f'{str_tmp}.prob.svg',
    )
    plot.save(
        fig,
        tup_size = tup_size,
        str_save = f'{str_tmp}.prob.pdf',
    )
#'''

plt.show()
