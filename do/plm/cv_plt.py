import numpy as np
import matplotlib.pyplot as plt

def def_plt(
    list_header: list,
    float_timescale: float
) -> None:

    int_nplot = len(list2d_header)
    fig, axs = plt.subplots(int_nplot, 1, sharex='all')
    if int_nplot==1:
        axs = [axs]

    with open('COLVAR', 'r') as colvar:
        list_header = colvar.readline().split()[2:]
    data = np.genfromtxt("COLVAR", dtype=None, names=list_header, invalid_raise=False)
    #print(data.dtype)

    for int_i in range(int_nplot):
        list_header = list2d_header[int_i]
        for str_header in list_header:
            if str_header not in dict_label:
                str_label = str_header
            else:
                str_label = dict_label[str_header]
            axs[int_i].scatter(data['time']*float_timescale, data[str_header], label=str_label, s=0.5)
            axs[int_i].legend(loc='upper left') 
    #axs[-1].set_xlabel('Time(ps)')
    axs[0].set_ylabel('CV')

def def_multi_plt(
    list_header: list,
    float_timescale: float
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
                str_header = list_header[int_i]
                if str_header not in data.dtype.names:
                    continue
                axs[int_i].scatter(data['time']*float_timescale, data[str_header], s=0.5)
            if str_line[0] != '#':
                break
            list_tmp = [str_line]
            list_field = str_line.split()[2:]
    for int_i in range(int_nplot):
        str_header = list_header[int_i]
        if str_header not in dict_label:
            str_label = str_header
        else:
            str_label = dict_label[str_header]
        axs[int_i].set_ylabel(str_label)

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
        'dist_vp_c',
        'dh2x_o_0_h',
        'metadbias',
        'dist_o_0_h',
        #'uw_dist_o_0_hbias',
        'cn_o_0_h',
        #'uw_cn_o_0_hbias',
        'cnn_o_w_hmin'
    ],
    float_timescale = 1/0.0005/40
)

plt.show()
