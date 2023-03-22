import pandas as pd
from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

def plt_state(
    float_xscale: float = 1,
    str_xlabel: str = None,
    tup_ylim: tuple = None,
    str_file: str = 'carbonic_rolling.csv',
):

    df_data = pd.read_csv(str_file)
    print(df_data)
    
    list_header = ['HCO3', 'TT', 'CT', 'CC']
    list_pos = [1, 2, 3, 4]
    
    fig, ax = plt.subplots()
    
    df_new = df_data.where(df_data.isnull(), 1)
    df_new['frame'] = df_data['frame']
    
    for idx, header in enumerate(list_header):
        df_data_tmp = df_data[df_data[header].notnull()]
        df_new_tmp = df_new[df_new[header].notnull()]
        if len(df_data_tmp)==0:
            continue
        ax.scatter( df_new_tmp['frame']*float_xscale, df_new_tmp[header]*list_pos[idx], s=2, edgecolors='none', c='tab:blue', alpha=df_data_tmp[header])

    ax.set_xlabel(str_xlabel)
    ax.set_yticks(list_pos)
    ax.set_yticklabels(['HCO$_3^-$', 'TT','CT','CC'])
    ax.set_ylim(tup_ylim)

    return fig, ax

'''
fig, ax = plt_state(
    float_xscale = 0.0004837769,
    str_xlabel = 'Time (ps)',
    tup_ylim = (0.5, 4.5),
)
plot.add_text(
    ax,
    dict_text = {'AIMD': (0.9, 0.9)}
)
plot.save(
    fig,
    tup_size = (8.6*cm, 3.5*cm),
    str_save = 'carbonic_rolling',
    list_type = ['png'],
)
#'''

#'''
fig, ax = plt_state(
    str_file = 'carbonic_rolling.csv',
)

plot.save(
    fig,
    tup_size = (8.6*cm, 3*cm),
)
#'''

'''
fig, ax = plt_state(
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
