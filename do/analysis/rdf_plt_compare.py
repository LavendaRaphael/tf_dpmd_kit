from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis

# data

dict_input = {
    'DPMD (400 ps)': ('/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/04.H2CO3_TT_H2O_126/330K/rdf/', '0050000_0130000'),
    'AIMD (28 ps)': ('/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/rdf/', '0001000_0057877'),
}

def gen_data(
    str_pair: str,
    dict_input: dict,
) -> dict:
    
    dict_data = {}
    for str_label, tup_data in dict_input.items():
        dict_data[str_label] = f'{tup_data[0]}/rdf.{str_pair}.{tup_data[1]}.csv'

    return dict_data

str_pair = 'c.o_w'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    tup_xrange = (2,6),
    tup_yrange = (0,2.5),
    str_save = f'rdf.{str_pair}.compare.pdf',
    str_title = r'C-O$_W$'
)

str_pair = 'o_1.h_w'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    tup_xrange = (1,6),
    tup_yrange = (0,2),
    str_save = f'rdf.{str_pair}.compare.pdf',
    str_title = r'$^=$O-H$_W$'
)

str_pair = 'o_0_2.h_w'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    tup_xrange = (1,6),
    tup_yrange = (0,2),
    str_save = f'rdf.{str_pair}.compare.pdf',
    str_title = r'O$_H$-H$_W$'
)

str_pair = 'h_0_1.o_w'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    tup_xrange = (1,6),
    tup_yrange = (0,3),
    str_save = f'rdf.{str_pair}.compare.pdf',
    str_title = r'H$_O$-O$_W$'
)

str_pair = 'o_w.o_w'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair=str_pair, dict_input=dict_input),
    tup_xrange = (2,6),
    tup_yrange = (0,3.5),
    str_save = f'rdf.{str_pair}.compare.pdf',
    str_title = r'O$_W$-O$_W$'
)

plt.show() 
