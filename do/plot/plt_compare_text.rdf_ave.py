from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

str_pair = 'o_1.h_w'
plot.plt_compare_text(
    dict_data = {
        '4-28ps': f'rdf.{str_pair}.ave.csv'
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    str_title = r'$^=$O-H$_W$',
    tup_xlim = (1,6),
    tup_ylim = (0,2.5),
    str_save = f'rdf.{str_pair}.ave.pdf',
    bool_error = True,
)
str_pair = 'h_0_1.o_w'
plot.plt_compare_text(
    dict_data = {
        '4-28ps': f'rdf.{str_pair}.ave.csv'
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    str_title = r'H$_T$-O$_W$',
    tup_xlim = (1,6),
    tup_ylim = (0,3),
    str_save = f'rdf.{str_pair}.ave.pdf',
    bool_error = True,
)
str_pair = 'o_0_2.h_w'
plot.plt_compare_text(
    dict_data = {
        '4-28ps': f'rdf.{str_pair}.ave.csv'
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    str_title = r'O$_T$-H$_W$',
    tup_xlim = (1,6),
    tup_ylim = (0,2.5),
    str_save = f'rdf.{str_pair}.ave.pdf',
    bool_error = True,
)
str_pair = 'o_w.o_w'
plot.plt_compare_text(
    dict_data = {
        '4-28ps': f'rdf.{str_pair}.ave.csv'
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    str_title = r'O$_W$-O$_W$',
    tup_xlim = (2,6),
    tup_ylim = (0,4),
    str_save = 'rdf.o_w.o_w.ave.pdf',
    bool_error = True,
)

plt.show()
