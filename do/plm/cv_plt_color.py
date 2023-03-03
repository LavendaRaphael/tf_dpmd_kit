from tf_dpmd_kit import plm
from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

dict_label = {
    'dist_vp_c': r'R(V$_{Px}$C)',
    'cn_o_h': r'CN(O$_{CA}$H)',
}

def gen_dict_header(
    list_header: list,
    dict_label: dict,
) -> dict:
    dict_header = {}
    for str_header in list_header:
        if str_header in dict_label:
            dict_header[str_header] = dict_label[str_header]
        else:
            dict_header[str_header] = str_header
    return dict_header

list_header = ['cn_o_h']

plot.set_rcparam()
cm = 1/2.54

#'''
fig, axs = plm.colvar_plt(
    list_data = ['COLVAR.0', 'COLVAR.1'],
    dict_header = gen_dict_header(list_header, dict_label),
    tup_xlim = (0, 2000),
    #tup_ylim = (0,15),
    str_xlabel = None,
    dict_color = {
        (     0, 400000): 'tab:blue',
    },
    str_color = 'tab:grey',
    float_scatters = 1.0,
)
ax = axs[0]
plot.save(
    fig,
    str_save = 'cn_o_h.time.color.png',
    tup_size = (8.6*cm, 3*cm),
)
#'''

plt.show()
