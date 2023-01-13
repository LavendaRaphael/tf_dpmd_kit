from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

dict_label = {
    'cn_o_h': r'CN(O$_C$H)',
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

#'''
fig, axs = plm.colvar_plt(
    dict_header = gen_dict_header(list_header, dict_label),
    tup_xlim = (0,None),
    tup_ylim = (0,None),
    str_label = 'AIMD'
)
ax = axs[0]
plm.insert_img(
    fig,ax,
    dict_img = {
        '10000.png': ( 0.0, 0, 0.2, 0.3),
        '26950.png': (0.25, 0, 0.2, 0.3),
        '60060.png': ( 0.5, 0, 0.2, 0.3),
        '63000.png': (0.75, 0, 0.2, 0.3),
    },
    dict_arrow = {
        (10000*0.0004837769, 1.7): ( 0.1, 0.3),
        (26950*0.0004837769, 1.0): (0.35, 0.3),
        (60060*0.0004837769, 1.2): ( 0.6, 0.3),
        (63000*0.0004837769, 1.7): (0.85, 0.3),
    },
    str_save = 'cn_o_h.time.img.png',
    tup_size = (11,5)
)
#'''


plt.show()
