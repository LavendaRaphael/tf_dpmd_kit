import matplotlib.pyplot as plt
from tf_dpmd_kit import plm

dict_label = {
}

#'''
list_header = [
    'dist_vp_o_1_2', 'dhx_o_0_h','dhx_o_2_h', 'dist_o_0_h','cn_o_0_h'
]
for str_tmp in list_header:
    if str_tmp in dict_label:
        str_label = dict_label[str_header]
    else:
        str_label = str_tmp
    plm.colvar_hist_plt(
        str_header = str_tmp,
        str_save = f'{str_tmp}.prob.pdf',
        str_label = str_tmp
    )
#'''
plt.show()
