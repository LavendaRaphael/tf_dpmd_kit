import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def def_multi_plt(
    list_header: list,
    float_timescale: float=1.0,
    tup_xlim: tuple = None,
    str_save: str = None,
    str_color: str = None,
) -> None:

    rc('font',**{'size':15, 'family':'sans-serif','sans-serif':['Arial']})

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
                axs[int_i].scatter(data['time']*float_timescale, data[str_header], s=0.5, color=str_color)
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
    axs[-1].set_xlabel(('time (ns)'))
    if str_save:
        fig.set_size_inches(11, 5)
        fig.savefig(str_save, bbox_inches='tight', dpi=300)

dict_label = {
    'dist_vp_o_1': r'R(V$_p$O$_1$)',
    'dist_vp_o_2': r'R(V$_p$O$_2$)',
    'dist_vp_o_1_2': r'R(V$_p$O$_C$)',
    'dist_vp_c': r'R(V$_p$C)',
    'cost_o_h': 'cos dh(OH)',
    'dist_o_0_h': r'R(O$_0$H)',
    'metadbias': 'V(s,t)',
    'metadrbias': 'M(s,t)', 
    'metadrct': 'c(t)', 
}

def_multi_plt(
    list_header = [
        #['dist_vp_o_1'],
        #['dist_vp_o_2'],
        ['dist_vp_o_1_2'],

        ['dist_vp_c'],
        #['dist_vp_c', (0.5, 1.5)],
        #['lw_vp_cbias'],

        ['cn_o_1_2_h'],

        #['dhx_o_0_h'],
        #['dhx_o_2_h'],

        ['metadbias'],
        #['metadrbias'],
        #['metadrct'],

        ['dist_o_0_h'],
        #['dist_o_0_h', (1.1, 1.2)],
        #['uw_dist_o_0_hbias'],

        ['cn_o_0_h'],
        #['uw_cn_o_0_hbias'],
    ],
    float_timescale = 1/1000,
    #float_timescale = 1/0.0005/50,
    #float_timescale = 1/0.000048378/10
    #tup_xlim = (360800, 363300),
    #str_save = 'cv.png'
)

plt.show()
