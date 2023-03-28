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
    legend_kw: dict = None,
    bool_legend: bool = True,
):

    if dict_label is None:
        dict_label = {}
    if dict_color is None:
        dict_color = {}
    if legend_kw is None:
        legend_kw = {}

    df_data = pd.read_csv(file_data)
    for header, marker in zip(list_header, Line2D.filled_markers):
        label = header
        if header in dict_label:
            label = dict_label[header]
        color = None
        if header in dict_color:
            color = dict_color[header]

        ax.errorbar(df_data['temperature(K)'], df_data[header], label=label, ls=':', marker=marker, markersize=2, lw=1, color=color)

    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(ylabel)
    if bool_legend:
        ax.legend(
            frameon = False,
            **legend_kw
        )

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
        file_data = 'carbonic_statistic.temperature_timemean.csv',
        list_header = ['TT', 'CT', 'CC','H2CO3', 'HCO3'],
        ylabel = 'Lifetime (ps)',
        dict_label = dict_label, 
        dict_color = dict_color,
        bool_legend = False,
    )
    plot.add_text(
        ax,
        dict_text = {
            (355, 1000): 'TT',
            (355,  500): r'H$_2$CO$_3$',
            (355,  100): 'CT',
            (355,   30): r'HCO$_3^-$',
            (355,   15): 'CC',
        }
    )
    ax.set_yscale('log')
    ax.set_xlim(None, 380)
    plot.save(
        fig,
        file_save = 'carbonic_statistic.temperature_timemean',
        list_type = ['pdf', 'svg']
    )

    fig, ax = plt.subplots(figsize=(4.3*cm, 4*cm))
    run(
        ax,
        file_data = 'carbonic_statistic.temperature_count.csv',
        list_header = ['TT', 'CT', 'CC', 'H2CO3'],
        ylabel = r'Count (ns$^{-1}$)',
        dict_label = dict_label,
        dict_color = dict_color,
    )
    plot.save(
        fig,
        file_save = 'carbonic_statistic.temperature_count',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
