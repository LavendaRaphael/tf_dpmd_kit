from matplotlib import pyplot as plt
from tf_dpmd_kit import plot
import numpy as np

def run(
    dict_data,
    str_pair,
    str_text,
    tup_ylim,
) -> None:

    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(4.3*cm, 3.7*cm))

    plot.plt_compare(
        ax,
        dict_data = dict_data, 
        float_lw = 0.5,
    )
    ax.set_xlabel('r (Å)')
    ax.set_ylabel('g(r)')
    ax.set_xlim(1,6)
    ax.set_ylim(tup_ylim)
    ax.legend(
        frameon = False,
        fontsize = 'small',
        labelspacing = 0.1
    )
    plot.add_text(
        ax,
        dict_text = {
            (0.2, 0.9): str_text,
        },
        transform = ax.transAxes
    )
    plot.save(
        fig,
        file_save = f'rdf.{str_pair}.converge',
        list_type = ['pdf', 'svg']
    )

def rdfcsv(file):

    data = np.loadtxt(file)
    return data[:,0], data[:,1]

def get_dict_data(pair):

    return {
        '0-0.5 ns'  : rdfcsv(f'rdf.{pair}.0000000_0100000.csv'),
        '0.5-1.0 ns': rdfcsv(f'rdf.{pair}.0100000_0200000.csv'),
        '1.0-1.5 ns': rdfcsv(f'rdf.{pair}.0200000_0300000.csv'),
        '1.5-2.0 ns': rdfcsv(f'rdf.{pair}.0300000_0400000.csv'),
    }

def main():

    plot.set_rcparam()

    str_pair = 'h_0_1.o_w'
    run(
        dict_data = get_dict_data(str_pair),
        str_text = r'H$_{OH}$-O$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 3.5),
    )
    str_pair = 'o_0_2.h_w'
    run(
        dict_data = get_dict_data(str_pair),
        str_text = r'O$_{OH}$-H$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 2.5),
    )
    str_pair = 'o_1.h_w'
    run(
        dict_data = get_dict_data(str_pair),
        str_text = r'$^=$O-H$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 2.5),
    )
    str_pair = 'o_w.o_w'
    run(
        dict_data = get_dict_data(str_pair),
        str_text = r'O$_W$-O$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 4),
    )
    
    plt.show()

main()
