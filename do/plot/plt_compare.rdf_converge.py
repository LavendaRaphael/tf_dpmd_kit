from matplotlib import pyplot as plt
from tf_dpmd_kit import plot

plot.set_rcparam()
cm = 1/2.54

def gen_label(
    float_timestep: float,
    tup_id: tuple,
):
    float_start = tup_id[0]*float_timestep
    float_end = tup_id[-1]*float_timestep
    return f'{float_start:.2f}-{float_end:.2f}ns'

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

def run(
    dict_title: dict,
    dict_data: dict,
    dict_xlim: dict,
    dict_ylim: dict,
    dict_save: dict,
) -> None:

    for str_key, str_titile in dict_title.items():
        fig, ax = plot.plt_compare(
            dict_data = dict_data[str_key],
            tup_xlim = dict_xlim[str_key],
            tup_ylim = dict_ylim[str_key],
            str_xlabel = 'r (Å)',
            str_ylabel = 'g(r)',
            float_lw = 1,
        )
        plot.add_text(
            ax,
            dict_text = { 
                dict_title[str_key]: (0.1, 0.9)
            }
        )
        plot.save(
            fig,
            tup_size = (4.3*cm, 4*cm),
            str_save = dict_save[str_key],
            list_type = ['pdf', 'svg']
        )

# setup

dict_title = {
    'o_0_2.h_w': r'O$_T$-H$_W$',
    'o_1.h_w': r'$^=$O-H$_W$',
    'h_0_1.o_w': r'H$_T$-O$_W$',
    'o_w.o_w': r'O$_W$-O$_W$',
}

dict_xlim = {
    'o_0_2.h_w': (1,6),
    'o_1.h_w': (1,6),
    'h_0_1.o_w': (1,6),
    'o_w.o_w': (2,6),
}

dict_ylim = {
    'o_0_2.h_w': (0,2.5),
    'o_1.h_w': (0,2.5),
    'h_0_1.o_w': (0,3),
    'o_w.o_w': (0,4.0)
}

list_id = [
    (      0,100000),
    ( 100000,200000),
    ( 200000,300000),
    ( 300000,400000),
]
float_timestep = 0.000005
dict_data = {
    'o_0_2.h_w': gen_data('o_0_2.h_w', list_id, float_timestep),
    'o_1.h_w': gen_data('o_1.h_w', list_id, float_timestep),
    'h_0_1.o_w': gen_data('h_0_1.o_w', list_id, float_timestep),
    'o_w.o_w': gen_data('o_w.o_w', list_id, float_timestep),
}

dict_save = {
    'o_0_2.h_w': f'rdf.o_0_2.h_w.converge',
    'o_1.h_w': f'rdf.o_1.h_w.converge',
    'h_0_1.o_w': f'rdf.h_0_1.o_w.converge',
    'o_w.o_w': f'rdf.o_w.o_w.converge',
}

run(
    dict_title = dict_title,
    dict_data = dict_data,
    dict_xlim = dict_xlim,
    dict_ylim = dict_ylim,
    dict_save = dict_save,
)

plt.show()

