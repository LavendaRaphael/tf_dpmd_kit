from tf_dpmd_kit import plm
from tf_dpmd_kit import analysis
import matplotlib.pyplot as plt

dict_color = plm.dict_color_temperature(
    tup_colormap = (280, 330),
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320,
        '330K': 330,
    }
)

str_pair = 'o_1.h_w'
analysis.plt_compare_text(
    dict_data = {
        '280K': f'280K/rdf/rdf.{str_pair}.ave.csv',
        '290K': f'290K/rdf/rdf.{str_pair}.ave.csv',
        '300K': f'300K/rdf/rdf.{str_pair}.ave.csv',
        '310K': f'310K/rdf/rdf.{str_pair}.ave.csv',
        '320K': f'320K/rdf/rdf.{str_pair}.ave.csv',
        '330K': f'330K/rdf/rdf.{str_pair}.ave.csv',
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    dict_color = dict_color,
    str_title = r'$^=$O-H$_W$',
    tup_xlim = (1,6),
    tup_ylim = (0,2.5),
    str_save = f'rdf.{str_pair}.ave.temperature.pdf',
    bool_error = True,
)
str_pair = 'h_0_1.o_w'
analysis.plt_compare_text(
    dict_data = {
        '280K': f'280K/rdf/rdf.{str_pair}.ave.csv',
        '290K': f'290K/rdf/rdf.{str_pair}.ave.csv',
        '300K': f'300K/rdf/rdf.{str_pair}.ave.csv',
        '310K': f'310K/rdf/rdf.{str_pair}.ave.csv',
        '320K': f'320K/rdf/rdf.{str_pair}.ave.csv',
        '330K': f'330K/rdf/rdf.{str_pair}.ave.csv',
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    dict_color = dict_color,
    str_title = r'H$_T$-O$_W$',
    tup_xlim = (1,6),
    tup_ylim = (0,3),
    str_save = f'rdf.{str_pair}.ave.temperature.pdf',
    bool_error = True,
)
analysis.plt_compare_text(
    dict_data = {
        '280K': '280K/rdf/rdf.o_w.o_w.ave.csv',
        '290K': '290K/rdf/rdf.o_w.o_w.ave.csv',
        '300K': '300K/rdf/rdf.o_w.o_w.ave.csv',
        '310K': '310K/rdf/rdf.o_w.o_w.ave.csv',
        '320K': '320K/rdf/rdf.o_w.o_w.ave.csv',
        '330K': '330K/rdf/rdf.o_w.o_w.ave.csv',
    },
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
    dict_color = dict_color,
    str_title = r'O$_W$-O$_W$',
    tup_xlim = (2,6),
    tup_ylim = (0,4),
    str_save = 'rdf.o_w.o_w.ave.temperature.pdf',
    bool_error = True,
)

plt.show()
