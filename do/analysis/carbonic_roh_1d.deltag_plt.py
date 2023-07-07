import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
from tf_dpmd_kit import plm
from tf_dpmd_kit import analysis
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.lines import Line2D

def run(
    ax,
):

    str_exp = '/home/faye/research/202203_MDCarbonicAcid/record/ref/'

    dict_data = {
        '2019 PNAS': str_exp+'/2019_PNAS_DanielAminov/Fig_1_deltag.csv',
        '2016 JPCB': str_exp+'/2016_JPCB_DinePines/deltag.csv',
        '2010 JPCA': str_exp+'/2010_JPCA_WangXiaoguang/Sfig_3_deltag.csv',
        '2009 Science': str_exp+'/2009_Science_KatrinAdamczyk/deltag.csv',
    }
    list_marker = Line2D.filled_markers[1:]
    marker_id = 0
    for label, file in dict_data.items():
        print(label)
        np_data = np.loadtxt(file, ndmin=2)
        ax.errorbar(np_data[:,0], np_data[:,1]/4.184, yerr=np_data[:,2]/4.184, ls='', label=label, capsize=2, lw=1, markersize=2, marker=list_marker[marker_id])
        marker_id += 1

    df = pd.read_csv('carbonic_roh_1d.pka.csv')
    df = df.groupby('state').get_group('H2CO3')
    ax.errorbar(df['temperature(K)'], df['deltag(kcal/mol)'], ls=':', marker=list_marker[marker_id], capsize=2, lw=1, markersize=2, label='this work')
    
    ax.set_xlabel(r'Temperature (K)')
    ax.set_ylabel('$\Delta$G (kcal/mol)')
    ax.legend()

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_roh_1d.deltag',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
