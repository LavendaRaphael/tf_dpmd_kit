from tf_dpmd_kit import plot
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

def run(
    dict_data,
    text,
    file_save,
    ylim,
):

    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(4.3*cm, 3.7*cm))
    for label, data in dict_data.items():
        ax.plot(data[0], data[1], label=label, lw=1)
    ax.text(
        x=0.1,
        y=0.9,
        s = text,
        transform=ax.transAxes
    )
    ax.set_xlim(1,6)
    ax.set_ylim(ylim)
    ax.set_xlabel('r (Å)')
    ax.legend()
    plot.save(
        fig,
        file_save = file_save,
        list_type = ['pdf', 'svg']
    )

def main():

    plot.set_rcparam()
    mpl.rcParams['figure.dpi'] = 300

    df_aimd = pd.read_csv('../../../../CPBO/CC/carbonicrdf/carbonicrdf.csv', index_col='r(ang)')
    df_dpmd = pd.read_csv('carbonicrdf.csv', index_col='r(ang)')

    pair = 'cc.o_nyl.h_w'
    run(
        dict_data = {
            'DPMD': [df_dpmd.index, df_dpmd[pair]],
            'AIMD': [df_aimd.index, df_aimd[pair]],
        },
        text = r'$^=$O-H$\mathregular{_W}$',
        file_save = f'rdf.{pair}',
        ylim = (0, 2)
    )
    pair = 'cc.o_oh.h_w'
    run(
        dict_data = {
            'DPMD': [df_dpmd.index, df_dpmd[pair]],
            'AIMD': [df_aimd.index, df_aimd[pair]],
        },
        text = r'O$\mathregular{_{OH}}$-H$\mathregular{_W}$',
        file_save = f'rdf.{pair}',
        ylim = (0, 2)
    )
    pair = 'cc.h_oh.o_w'
    run(
        dict_data = {
            'DPMD': [df_dpmd.index, df_dpmd[pair]],
            'AIMD': [df_aimd.index, df_aimd[pair]],
        },
        text = r'H$\mathregular{_{OH}}$-O$\mathregular{_W}$',
        file_save = f'rdf.{pair}',
        ylim = (0, 3)
    )
    pair = 'o_w.o_w'
    run(
        dict_data = {
            'DPMD': [df_dpmd.index, df_dpmd[pair]],
            'AIMD': [df_aimd.index, df_aimd[pair]],
        },
        text = r'O$\mathregular{_W}$-O$\mathregular{_W}$',
        file_save = f'rdf.{pair}',
        ylim = (0, 4)
    )
    plt.show()

main()
