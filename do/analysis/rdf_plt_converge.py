from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis

def gen_label(
    float_timestep: float,
    tup_id: tuple,
):
    int_start = int(tup_id[0]*float_timestep)
    int_end = int(tup_id[-1]*float_timestep)
    return f'{int_start}-{int_end}ps'

def gen_data(
    str_pair: str,
    list_id: list,
    float_timestep: float
) -> dict:

    dict_data = {}
    for tup_id in list_id:
        str_label = gen_label(float_timestep, tup_id)
        dict_data[str_label] = f'rdf.{str_pair}.{tup_id[0]:07d}_{tup_id[-1]:07d}.csv'

    return dict_data

# setup

#'''
list_id = [
    (      0, 50000),
    (  50000,100000),
    ( 100000,150000),
    ( 150000,200000)
]
float_timestep = 0.005
#'''

str_pair = 'o_1.h_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair, list_id, float_timestep),
    tup_xlim = (1,6),
    tup_ylim = (0,2),
    str_save = f'rdf.{str_pair}.converge.pdf',
    str_title = r'$^=$O-H$_W$',
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
)
str_pair = 'o_0_2.h_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair, list_id, float_timestep),
    tup_xlim = (1,6),
    tup_ylim = (0,2),
    str_save = f'rdf.{str_pair}.converge.pdf',
    str_title = r'O$_H$-H$_W$',
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
)
str_pair = 'h_0_1.o_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair, list_id, float_timestep),
    tup_xlim = (1,6),
    tup_ylim = (0,3),
    str_save = f'rdf.{str_pair}.converge.pdf',
    str_title = r'H$_O$-O$_W$',
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
)
str_pair = 'o_w.o_w.TT'
analysis.rdf_plt_compare(
    dict_data = gen_data(str_pair, list_id, float_timestep),
    tup_xlim = (2,6),
    tup_ylim = (0,4),
    str_save = f'rdf.{str_pair}.converge.pdf',
    str_title = r'O$_W$-O$_W$',
    str_xlabel = 'r (Å)',
    str_ylabel = 'g(r)',
)
plt.show()

