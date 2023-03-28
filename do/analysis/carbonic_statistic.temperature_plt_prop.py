import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from tf_dpmd_kit import plot
import pandas as pd


def run(
    ax,
    file_data: str,
    list_header: list,
    ylabel: str,
    dict_label: dict = None,
    dict_color: dict = None,
    ylim: tuple = None,
):

    if dict_label is None:
        dict_label = {}
    if dict_color is None:
        dict_color = {}

    df_data = pd.read_csv(file_data)
    ser_0 = pd.Series([0]*len(df_data))
    for header, marker in zip(list_header, Line2D.filled_markers):
        label = header
        if header in dict_label:
            label = dict_label[header]
        color = None
        if header in dict_color:
            color = dict_color[header]
        ser_1 = ser_0 + df_data[header]
        ax.fill_between(df_data['temperature(K)'], ser_0, ser_1, label=label, lw=1, color=color, alpha=0.5)
        ser_0 = ser_1

    plot.add_text(
        ax,
        dict_text = {
            (330, 0.75): 'TT',
            (330, 0.85): 'CT',
            (330, 0.90): 'CC',
            (330, 0.97): r'HCO$_3^-$',
        },
        va = 'center',
        ha = 'center',
        color = 'white',
        fontweight = 'bold',
        fontsize = 8,
    ),
    plot.add_arrow(
        ax,
        dict_arrow = {
            'x': [(335, 0.90), (347, 0.92)]
        },
        arrowstyle = 'simple, head_length=2, head_width=2, tail_width=0.1',
        color = 'white'
    )

    ax.set_ylim(ylim)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(ylabel)

def main():

    plot.set_rcparam()
    cm = 1/2.54

    dict_label = {
        'H2CO3': r'H$_2$CO$_3$',
        'HCO3': r'HCO$_3^-$',
    }
    dict_color = {
        'TT': 'tab:blue',
        'CT': 'tab:orange',
        'CC': 'tab:green',
        'H2CO3': 'tab:red',
        'HCO3': 'tab:purple',
    }

    fig, ax = plt.subplots(figsize=(4.3*cm, 4*cm))
    run(
        ax,
        file_data = 'carbonic_statistic.temperature_prop.csv',
        list_header = ['TT', 'CT', 'CC', 'HCO3'],
        ylabel = 'Proportion',
        dict_label = dict_label,
        dict_color = dict_color,
        ylim = (0.7, 1),
    )
    plot.save(
        fig,
        file_save = 'carbonic_statistic.temperature_prop',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
