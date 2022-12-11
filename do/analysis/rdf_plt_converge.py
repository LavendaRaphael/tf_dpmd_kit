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

def plot(
    dict_title: dict,
    dict_data: dict,
    dict_xlim: dict,
    dict_ylim: dict,
    dict_save: dict,
) -> None:

    for str_key, str_titile in dict_title.items():
        analysis.rdf_plt_compare(
            dict_data = dict_data[str_key],
            tup_xlim = dict_xlim[str_key],
            tup_ylim = dict_ylim[str_key],
            str_save = dict_save[str_key],
            str_title = dict_title[str_key],
            str_xlabel = 'r (Å)',
            str_ylabel = 'g(r)',
        )

# setup

dict_title = {
    'o_0.h_w': r'O$_H$-H$_W$',
    'o_1_2.h_w': r'$^=$O-H$_W$',
    'h_0.o_w': r'H$_O$-O$_W$',
}

dict_xlim = {
    'o_0.h_w': (1,6),
    'o_1_2.h_w': (1,6),
    'h_0.o_w': (1,6),
}

dict_ylim = {
    'o_0.h_w': (0,2.5),
    'o_1_2.h_w': (0,2.5),
    'h_0.o_w': (0,3),
}

list_id = [
    (      0,100000),
    ( 100000,200000),
    ( 200000,300000),
    ( 300000,400000)
]
float_timestep = 0.005
dict_data = {
    'o_0.h_w': gen_data('o_0.h_w', list_id, float_timestep),
    'o_1_2.h_w': gen_data('o_1_2.h_w', list_id, float_timestep),
    'h_0.o_w': gen_data('h_0.o_w', list_id, float_timestep),
}

dict_save = {
    'o_0.h_w': f'rdf.o_0.h_w.converge.pdf',
    'o_1_2.h_w': f'rdf.o_1_2.h_w.converge.pdf',
    'h_0.o_w': f'rdf.h_0.o_w.converge.pdf',
}

plot(
    dict_title = dict_title,
    dict_data = dict_data,
    dict_xlim = dict_xlim,
    dict_ylim = dict_ylim,
    dict_save = dict_save,
)

plt.show()

