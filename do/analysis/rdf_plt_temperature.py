from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis,plm

# data

dict_input = {
    '280K': ('280K.0/rdf/', '0050000_0200000'),
    '290K': ('290K.0/rdf/', '0050000_0200000'),
    '300K': ('300K.0/rdf/', '0050000_0200000'),
    '310K': ('310K.0/rdf/', '0050000_0200000'),
    '320K': ('320K.0/rdf/', '0050000_0200000'),
}

dict_temperature = {
    '280K': 280,
    '290K': 290,
    '300K': 300,
    '310K': 310,
    '320K': 320
}

tup_colormap = (280,330)

def gen_data(
    str_pair: str,
    dict_input: dict,
) -> dict:
    
    dict_data = {}
    for str_label, tup_data in dict_input.items():
        dict_data[str_label] = f'{tup_data[0]}/rdf.{str_pair}.{tup_data[1]}.csv'

    return dict_data

str_pair = 'o_1.h_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    dict_temperature = dict_temperature,
    tup_xlim = (1,6),
    tup_ylim = (0,2.5),
    str_save = f'rdf.{str_pair}.pdf',
    str_title = r'$^=$O-H$_W$',
    tup_colormap = tup_colormap
)

str_pair = 'o_0_2.h_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    dict_temperature = dict_temperature,
    tup_xlim = (1,6),
    tup_ylim = (0,2.5),
    str_save = f'rdf.{str_pair}.pdf',
    str_title = r'O$_H$-H$_W$',
    tup_colormap = tup_colormap
)

str_pair = 'h_0_1.o_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    dict_temperature = dict_temperature,
    tup_xlim = (1,6),
    tup_ylim = (0,3),
    str_save = f'rdf.{str_pair}.pdf',
    str_title = r'H$_O$-O$_W$',
    tup_colormap = tup_colormap
)

str_pair = 'o_w.o_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    dict_temperature = dict_temperature,
    tup_xlim = (2,6),
    tup_ylim = (0,4),
    str_save = f'rdf.{str_pair}.pdf',
    str_title = r'O$_W$-O$_W$',
    tup_colormap = tup_colormap
)

plt.show() 
