import pandas as pd
from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

def plt_state(
    float_xscale: float = 1,
    str_xlabel: str = None,
    tup_ylim: tuple = None,
):

    df_data = pd.read_csv('carbonic_rolling.csv')
    print(df_data)
    
    list_ticks = [0, 0.5, 1, 1.5, 1.8, 2, 2.2, 2.5, 3]
    
    fig, ax = plt.subplots()
    
    df_new = df_data.where(df_data.isnull(), 1)
    df_new['frame'] = df_data['frame']
    
    for idx, header in enumerate(df_data.columns[1:]):
        df_data_tmp = df_data[df_data[header].notnull()]
        df_new_tmp = df_new[df_new[header].notnull()]
        if len(df_data_tmp)==0:
            continue
        ax.scatter( df_new_tmp['frame']*float_xscale, df_new_tmp[header]*list_ticks[idx], s=2, edgecolors='none', c='tab:blue', alpha=df_data_tmp[header])

    ax.set_xlabel(str_xlabel)
    ax.set_yticks([ 1, 1.5, 1.8, 2, 2.2])
    ax.set_yticklabels(['HCO$_3^-$', '', 'TT','CT','CC'])
    ax.set_ylim(tup_ylim)

    return fig, ax

#'''
fig, ax = plt_state(
    float_xscale = 0.0004837769,
    str_xlabel = 'Time (ps)',
    tup_ylim = (0.7, 2.5),
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

'''
fig, ax = plt_state(
    float_xscale = 0.000005,
)

plot.save(
    fig,
    tup_size = (8.6*cm, 3*cm),
    str_save = 'carbonic_rolling',
    list_type = ['png'],
)
#'''

'''
fig, ax = plt_state(
    float_xscale = 0.000005,
    str_xlabel = 'Time (ns)',
    tup_ylim = (0.7, 2.5),
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
