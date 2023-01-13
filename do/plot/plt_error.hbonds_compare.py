from tf_dpmd_kit import plot
import matplotlib.pyplot as plt
import os

str_home = os.environ['homedir']

plot.plt_error(
    dict_data = {
        'TT': '../H2CO3_TT_H2O_126.0/hbonds.h_0_1.o_w.ave.temperature.scale.csv',
        'CT-T': 'hbonds.h_1.o_w.ave.temperature.csv',
        'CT-C': 'hbonds.h_0.o_w.ave.temperature.csv',
    },
    str_save = f'hbonds.ave.temperature.compare.pdf',
    str_xlabel = r'Temperature (K)',
    str_ylabel = 'n(H-bonds)'
)

plt.show()
