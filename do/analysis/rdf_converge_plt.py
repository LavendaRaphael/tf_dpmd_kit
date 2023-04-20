from matplotlib import pyplot as plt
from tf_dpmd_kit import plot

plot.set_rcparam()
cm = 1/2.54

def run(
    dict_data,
    str_pair,
    str_text,
) -> None:

    fig, ax = plot.plt_compare(
        dict_data = dict_data, 
        str_xlabel = 'r (Å)',
        str_ylabel = 'g(r)',
        float_lw = 1,
        legend_title = str_text,
        tup_xlim = (1,6),
    )
    plot.save(
        fig,
        tup_size = (4.3*cm, 3.7*cm),
        str_save = f'rdf.{str_pair}.converge',
        list_type = ['pdf', 'svg']
    )

str_pair = 'cc_h_oh.o_w'
run(
    dict_data = {
        '4-6ns':  f'rdf.{str_pair}.0800000_1200000.csv',
        '6-8ns':  f'rdf.{str_pair}.1200000_1600000.csv',
        '8-10ns': f'rdf.{str_pair}.1600000_2000000.csv',
    },
    str_text = r'H$_{OH}$-O$_W$',
    str_pair = str_pair,
)
str_pair = 'cc_o_oh.h_w'
run(
    dict_data = {
        '4-6ns':  f'rdf.{str_pair}.0800000_1200000.csv',
        '6-8ns':  f'rdf.{str_pair}.1200000_1600000.csv',
        '8-10ns': f'rdf.{str_pair}.1600000_2000000.csv',
    },
    str_text = r'O$_{OH}$-H$_W$',
    str_pair = str_pair,
)
str_pair = 'cc_o_c.h_w'
run(
    dict_data = {
        '4-6ns':  f'rdf.{str_pair}.0800000_1200000.csv',
        '6-8ns':  f'rdf.{str_pair}.1200000_1600000.csv',
        '8-10ns': f'rdf.{str_pair}.1600000_2000000.csv',
    },
    str_text = r'$^=$O-H$_W$',
    str_pair = str_pair,
)
str_pair = 'o_w.o_w'
run(
    dict_data = {
        '4-6ns':  f'rdf.{str_pair}.0800000_1200000.csv',
        '6-8ns':  f'rdf.{str_pair}.1200000_1600000.csv',
        '8-10ns': f'rdf.{str_pair}.1600000_2000000.csv',
    },
    str_text = r'O$_W$-O$_W$',
    str_pair = str_pair,
)

plt.show()
