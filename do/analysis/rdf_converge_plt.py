from matplotlib import pyplot as plt
from tf_dpmd_kit import plot

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
        str_xlabel = 'r (Å)',
        str_ylabel = 'g(r)',
        float_lw = 0.5,
        tup_xlim = (1,6),
        tup_ylim = tup_ylim,
        legend_kw = {
            'fontsize': 'small',
            'labelspacing': 0.1
        }
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

def main():

    plot.set_rcparam()

    str_pair = 'cc_h_oh.o_w'
    run(
        dict_data = {
            '0-6 ps'   : f'rdf.{str_pair}.0000000_0002497.csv',
            '6-12 ns'  : f'rdf.{str_pair}.0002497_0004978.csv',
            '12-18 ps' : f'rdf.{str_pair}.0004978_0007458.csv',
            '18-24 ps' : f'rdf.{str_pair}.0007458_0009939.csv',
            '24-30 ps' : f'rdf.{str_pair}.0009939_0012419.csv',
        },
        str_text = r'H$_{OH}$-O$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 3.5),
    )
    str_pair = 'cc_o_oh.h_w'
    run(
        dict_data = {
            '0-6 ps'   : f'rdf.{str_pair}.0000000_0002497.csv',
            '6-12 ns'  : f'rdf.{str_pair}.0002497_0004978.csv',
            '12-18 ps' : f'rdf.{str_pair}.0004978_0007458.csv',
            '18-24 ps' : f'rdf.{str_pair}.0007458_0009939.csv',
            '24-30 ps' : f'rdf.{str_pair}.0009939_0012419.csv',
        },
        str_text = r'O$_{OH}$-H$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 2.5),
    )
    str_pair = 'cc_o_c.h_w'
    run(
        dict_data = {
            '0-6 ps'   : f'rdf.{str_pair}.0000000_0002497.csv',
            '6-12 ns'  : f'rdf.{str_pair}.0002497_0004978.csv',
            '12-18 ps' : f'rdf.{str_pair}.0004978_0007458.csv',
            '18-24 ps' : f'rdf.{str_pair}.0007458_0009939.csv',
            '24-30 ps' : f'rdf.{str_pair}.0009939_0012419.csv',
        },
        str_text = r'$^=$O-H$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 2.5),
    )
    str_pair = 'o_w.o_w'
    run(
        dict_data = {
            '0-6 ps'   : f'rdf.{str_pair}.0000000_0002497.csv',
            '6-12 ns'  : f'rdf.{str_pair}.0002497_0004978.csv',
            '12-18 ps' : f'rdf.{str_pair}.0004978_0007458.csv',
            '18-24 ps' : f'rdf.{str_pair}.0007458_0009939.csv',
            '24-30 ps' : f'rdf.{str_pair}.0009939_0012419.csv',
        },
        str_text = r'O$_W$-O$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 4),
    )
    
    plt.show()

main()
