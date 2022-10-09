import numpy as np
import matplotlib.pyplot as plt

def def_multi_plt(
    list_header: list,
    float_timescale: float,
    tup_xlim: tuple = None,
    str_save: str = None,
) -> None:

    int_nplot = len(list_header)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as file_open:
        str_line = file_open.readline()
        list_tmp = [str_line]
        list_field = str_line.split()[2:]
        int_count = 0
        for str_line in file_open:
            print(list_tmp)
            list_tmp.append(str_line)
            for str_line in file_open:
                if str_line[:9] != '#! FIELDS':
                    list_tmp.append(str_line)
                else:
                    break
            data = np.genfromtxt(list_tmp, dtype=None, names=list_field, invalid_raise=False)
            for int_i in range(int_nplot):
                str_header = list_header[int_i][0]
                if str_header not in data.dtype.names:
                    continue
                axs[int_i].scatter(data['time']*float_timescale, data[str_header], s=0.5)
            if str_line[0] != '#':
                break
            list_tmp = [str_line]
            list_field = str_line.split()[2:]
    for int_i in range(int_nplot):
        str_header = list_header[int_i][0]
        if str_header not in dict_label:
            str_label = str_header
        else:
            str_label = dict_label[str_header]
        axs[int_i].set_ylabel(str_label)
        axs[int_i].set_xlim(tup_xlim)
        if len(list_header[int_i]) > 1:
            axs[int_i].set_ylim(list_header[int_i][1])
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

dict_label = {
    'dist_vp_o_1': 'R(O\N{SUBSCRIPT ONE}V\N{Latin Subscript Small Letter P})',
    'dist_vp_o_2': 'R(O\N{SUBSCRIPT TWO}V\N{Latin Subscript Small Letter P})',
    'dist_vp_c': 'R(CV\N{Latin Subscript Small Letter P})',
    'cost_o_h': 'cos dh(OH)',
    'dist_o_0_h': 'R(O\N{SUBSCRIPT ZERO}H)',
    'cost_o_0_h': 'cos dh(O\N{SUBSCRIPT ZERO}H)', 
}

def_multi_plt(
    list_header = [
        ['dist_vp_c', (0.5, 1.5)],
        #('dh2x_o_0_h'),
        ['metadbias'],
        #('dist_o_0_h'),
        #('uw_dist_o_0_hbias'),
        #('cn_o_0_h'),
        #('uw_cn_o_0_hbias'),
        #('cnn_o_w_hmin')
    ],
    float_timescale = 1/0.0005/40,
    tup_xlim = (360800, 363300),
    str_save = 'cv.pdf'
)

plt.show()
