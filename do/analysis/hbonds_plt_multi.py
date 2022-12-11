from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis

analysis.plt_subplots(
    dict_title = {
        'h_0_1.o_w': r'H$_T$-O$_W$',
        'h_w.o_0_2': r'O$_T$-H$_W$',
        'h_w.o_1': r'$^=$O-H$_W$',
        'h_w.o_w': r'O$_W$-H$_W$',
    },
    dict2d_data = {
        r'H$_2$CO$_3$-TT': {
            'h_0_1.o_w': 'hbonds.h_0_1.o_w.0000000_0057877.csv',
            'h_w.o_0_2': 'hbonds.h_w.o_0_2.0000000_0057877.csv',
            'h_w.o_1': 'hbonds.h_w.o_1.0000000_0057877.csv',
            'h_w.o_w': 'hbonds.h_w.o_w.0000000_0057877.csv'
        }
    },
    str_xlabel = 'Time (ps)',
    str_ylabel = 'N(HB)',
    str_save = f'hbonds.time.pdf',
)

plt.show()
