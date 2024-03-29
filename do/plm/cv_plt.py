from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

dict_label = {
    'dist_vp_o_1': r'R(V$_P$O$_1$)',
    'dist_vp_o_2': r'R(V$_P$O$_2$)',
    'dist_vp_o_1_2': r'R(V$_P$O$_C$)',
    'dist_vp_c': r'R(V$_P$C)',

    'dist_o_0_h': r'R(O$_0$H$_0$)',
    'dist_o_2_h': r'R(O$_2$H$_1$)',

    'cn_o_h': r'CN(O$_C$H)',
    'cn_o_0_h': r'CN(O$_0$H$_w$)',

    'dhx_o_0_h': 'dhx(O$_0$H$_0$)',
    'dhx_o_2_h': 'dhx(O$_2$H$_1$)',

    'metadbias': 'V(s,t)',
    'metadrbias': 'M(s,t)', 
    'metadrct': 'c(t)', 
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

list_header = [
    'cn_o_h',
    #'dist_o_0_h',
    #'dist_o_1_h',
    #'dhx_o_0_h',
    #'dhx_o_1_h',
    #'cos_o_0_h',
    #'cos_o_1_h',
    #'dhx_o_h',
]
plm.colvar_plt(
    list_data = ['COLVAR.0', 'COLVAR'],
    dict_header = gen_dict_header(list_header, dict_label),
    str_save = 'cv.time.png',
)

plt.show()
