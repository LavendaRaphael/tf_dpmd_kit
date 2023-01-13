from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis

def gen_label(
    float_timestep: float,
    tup_id: tuple,
):
    float_start = tup_id[0]*float_timestep
    float_end = tup_id[-1]*float_timestep
    return f'{float_start:.1f}-{float_end:.1f}ns'

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
    dict2d_data: dict,
    dict_xlim: dict,
    dict_ylim: dict,
    dict_save: dict,
) -> None:

    for str_key, str_titile in dict_title.items():
        analysis.plt_compare_text(
            dict_data = dict2d_data[str_key],
            tup_xlim = dict_xlim[str_key],
            tup_ylim = dict_ylim[str_key],
            str_save = dict_save[str_key],
            str_title = dict_title[str_key],
            str_xlabel = 'r (Å)',
            str_ylabel = 'g(r)',
        )

def gen_dict2d_data(
    list_pair: list,
    list_id: list,
    float_timestep: float,
) -> dict:

    dict2d_data = {}

    for str_pair in list_pair:
        dict2d_data[str_pair] = gen_data(
            str_pair = str_pair,
            list_id = list_id,
            float_timestep = float_timestep,
        )

    return dict2d_data

def gen_dict_save(
    list_pair: list,
) -> dict:

    dict_save = {}

    for str_pair in list_pair:
        dict_save[str_pair] = f'rdf.{str_pair}.converge.pdf'

    return dict_save

# setup

dict_title = {
    'h_0.o_w': r'H$_C$-O$_W$',
    'h_1.o_w': r'H$_T$-O$_W$',
    'o_0.h_w': r'O$_C$-H$_W$',
    'o_1.h_w': r'O$_T$-H$_W$',
    'o_2.h_w': r'$^=$O-H$_W$',
    'o_w.o_w': r'O$_W$-O$_W$'
}

dict_xlim = {
    'h_0.o_w': (1,6),
    'h_1.o_w': (1,6),
    'o_0.h_w': (1,6),
    'o_1.h_w': (1,6),
    'o_2.h_w': (1,6),
    'o_w.o_w': (2,6),
}

dict_ylim = {
    'h_0.o_w': (0,3.5),
    'h_1.o_w': (0,3.5),
    'o_0.h_w': (0,2.5),
    'o_1.h_w': (0,2.5),
    'o_2.h_w': (0,2.5),
    'o_w.o_w': (0,4.0),
}

list_pair = [
    'h_0.o_w',
    'h_1.o_w',
    'o_0.h_w',
    'o_1.h_w',
    'o_2.h_w',
    'o_w.o_w',
]
list_id = [
    (      0,100000),
    ( 100000,200000),
    ( 200000,300000),
    ( 300000,400000)
]
float_timestep = 0.000005
dict2d_data = gen_dict2d_data(
    float_timestep = float_timestep,
    list_pair = list_pair,
    list_id = list_id,
)

dict_save = gen_dict_save(
    list_pair = list_pair
)

plot(
    dict_title = dict_title,
    dict2d_data = dict2d_data,
    dict_xlim = dict_xlim,
    dict_ylim = dict_ylim,
    dict_save = dict_save,
)

plt.show()

