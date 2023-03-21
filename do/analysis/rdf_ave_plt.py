from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

plot.set_rcparam()
cm = 1/2.54

str_pair = 'o_0_1_2.h'
fig, ax = plot.plt_compare(
    dict_data = {
        'x': f'rdf.{str_pair}.ave.csv'
    },
    tup_xlim = (0.8, 2.6),
    tup_ylim = (0, 1),
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    bool_error = True,
    bool_legend = False,
    float_lw = 1,
)
plot.add_text(
    ax,
    dict_text = {
        r'O$_{CA}$-H': (0.4, 0.9)
    }
)
plot.save(
    fig,
    tup_size = (8.6*cm, 5*cm),
    str_save = f'rdf.{str_pair}.ave',
    list_type = ['pdf', 'svg'],
)

plt.show()
