from matplotlib import pyplot as plt
from tf_dpmd_kit import plot

plot.set_rcparam()
cm = 1/2.54

def run(
) -> None:

    fig, ax = plot.plt_compare(
        dict_data = {
            '5-10ns': 'rdf.o_0_1_2.h.1000000_2000000.csv',
            '10-15ns': 'rdf.o_0_1_2.h.2000000_3000000.csv',
            '15-20ns': 'rdf.o_0_1_2.h.3000000_4000000.csv',
            '20-25ns': 'rdf.o_0_1_2.h.4000000_5000000.csv',
        },
        str_xlabel = 'r (Å)',
        str_ylabel = 'g(r)',
        float_lw = 1,
    )
    plot.add_text(
        ax,
        dict_text = { 
            r'O$_{CA}$-H': (0.1, 0.9)
        }
    )
    plot.save(
        fig,
        tup_size = (8.6*cm, 4*cm),
        str_save = 'rdf.o_0_1_2.h.converge',
        list_type = ['pdf', 'svg']
    )

    plt.show()

run()


