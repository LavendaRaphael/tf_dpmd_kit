import pandas as pd
from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

def plt_state(
    float_xscale,
    str_xlabel,
    tup_ylim: tuple = None,
):

    df = pd.read_csv('carbonic_state.csv')
    print(df)
    
    list_ticks = [0, 0.5, 1, 1.5, 1.8, 2, 2.2, 2.5, 3]
    
    fig, ax = plt.subplots()
    for idx, header in enumerate(df.columns[1:]):
        ax.scatter( df['frame']*float_xscale, df[header]*list_ticks[idx], s=2, edgecolors='none', color='tab:blue')
    
    ax.set_xlabel(str_xlabel)
    ax.set_yticks(list_ticks)
    ax.set_yticklabels([r'CO$_3^{2-}$', '', 'HCO$_3^-$', '', 'TT','CT','CC', '', 'H$_3$CO$_3$'])
    ax.set_ylim(tup_ylim)

    return fig, ax

'''
fig, ax = plt_state(
    float_xscale = 0.0004837769,
    str_xlabel = 'Time (ps)',
)
'''

fig, ax = plt_state(
    float_xscale = 0.000005,
    str_xlabel = 'Time (ns)',
)

plot.save(
    fig,
    tup_size = (8.6*cm, 3.5*cm),
    str_save = 'carbonic',
    list_type = ['png'],
)

'''
df = pd.read_csv('carbonic.csv')
print(df)
fig, ax = plt.subplots()
ax.scatter( df['frame']*float_xscale, df['dihedral0(rad)'], s=0.5, edgecolors='none')
ax.scatter( df['frame']*float_xscale, df['dihedral1(rad)'], s=0.5, edgecolors='none')
ax.set_xlabel('Time (ps)')
plot.save(
    fig,
    tup_size = (8.6*cm, 3*cm),
    str_save = 'dihedral',
    list_type = ['png'],
)
#'''

plt.show()
