from tf_dpmd_kit import plot
from tf_dpmd_kit import analysis
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

fig, ax = plt.subplots(figsize=(8.6*cm, 3.5*cm))

#'''
analysis.carbonic_rolling_plt(
    ax,
    float_xscale = 0.0004837769,
    str_xlabel = 'Time (ps)',
    file_data = 'carbonic_state.csv',
    int_window = 100,
    list_header = ['CO3','HCO3','CC','CT','TT','H3CO3'],
    list_ypos = [0, 1, 2, 3, 4, 5],
    list_yticklabels = [r'CO$_3^{2-}$', 'HCO$_3^-$', 'CC','CT','TT', r'H$_3$CO$_3$'],
)


plot.save(
    fig,
    file_save = 'carbonic_rolling',
    list_type = ['pdf', 'svg']
)
#'''

'''
analysis.carbonic_rolling_plt(
    ax,
    float_xscale = 0.00001,
    str_xlabel = 'Time (ns)',
    file_data = 'carbonic_state.csv',
    int_window = 100,
    list_header = ['CO3','HCO3','CC','CT','TT','H3CO3'],
    list_ypos = [0, 1, 2, 3, 4, 5],
    list_yticklabels = [r'CO$_3^{2-}$', 'HCO$_3^-$', 'CC','CT','TT', r'H$_3$CO$_3$'],
)

#'''

plt.show()
