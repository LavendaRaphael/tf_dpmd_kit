from tf_dpmd_kit import plot
from tf_dpmd_kit import analysis
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

fig, ax = plt.subplots(figsize=(8.6*cm, 3.5*cm))

'''
analysis.carbonic_rolling_plt(
    float_xscale = 0.0004837769,
    str_xlabel = 'Time (ps)',
    tup_ylim = (0.5, 4.5),
    list_header = ['HCO3', 'TT', 'CT', 'CC'],
    list_ypos = [1, 2, 3, 4],
    list_yticklabels = ['HCO$_3^-$', 'TT','CT','CC'],
)
plot.add_text(
    ax,
    dict_text = {'AIMD': (0.9, 0.9)}
)
plot.save(
    fig,
    str_save = 'carbonic_rolling',
    list_size = ['pdf', 'svg']
)
#'''

#'''
analysis.carbonic_rolling_plt(
    ax,
    float_xscale = 0.00001,
    str_xlabel = 'Time (ns)',
    file_data = 'carbonic_state.csv',
    int_window = 100,
    list_header = ['CO3','HCO3','TT','CT','CC','H3CO3'],
    list_ypos = [0, 1, 2, 3, 4, 5],
    list_yticklabels = [r'CO$_3^{2-}$', 'HCO$_3^-$', 'TT','CT','CC', r'H$_3$CO$_3$'],
)

#'''

'''
analysis.carbonic_rolling_plt(
    ax,
    float_xscale = 0.000005,
    str_xlabel = 'Time (ns)',
    tup_ylim = (0.5, 4.5),
)
plot.add_text(
    ax,
    dict_text = {'DPMD': (0.9, 0.9)}
)
plot.save(
    fig,
    tup_size = (8.6*cm, 3.5*cm),
    str_save = 'carbonic_rolling',
    list_type = ['png'],
)
#'''
plt.show()
