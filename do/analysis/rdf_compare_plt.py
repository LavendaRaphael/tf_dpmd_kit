from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

def run(
    dict_data,
    str_pair,
    str_text,
    tup_ylim,
) -> None:

    cm = 1/2.54
    fig, ax = plt.subplots(figsize=(4.3*cm, 3.7*cm))

    plot.plt_compare(
        ax,
        dict_data = dict_data,
        str_xlabel = 'r (Å)',
        str_ylabel = 'g(r)',
        float_lw = 0.5,
        tup_xlim = (1,6),
        tup_ylim = tup_ylim,
        legend_kw = {
            'fontsize': 'small',
            'labelspacing': 0.1
        },
        bool_error = True,
    )
    plot.add_text(
        ax,
        dict_text = {
            (0.2, 0.9): str_text,
        },
        transform = ax.transAxes
    )
    plot.save(
        fig,
        file_save = f'rdf.{str_pair}.compare',
        list_type = ['pdf', 'svg']
    )

def main():

    plot.set_rcparam()

    dir_aimd = '/home/faye/research_d/202203_MDCarbonicAcid/server/07.md_water62/CPBO/CC/rdf/'
    dir_dpmd = '/home/faye/research_d/202203_MDCarbonicAcid/server/07.md_water62/DPMD/330K/CC/rdf/'

    str_pair = 'cc_h_oh.o_w'
    run(
        dict_data = {
            'AIMD': f'{dir_aimd}rdf.{str_pair}.ave.csv',
            'DPMD': f'{dir_dpmd}rdf.{str_pair}.ave.csv',
        },
        str_text = r'H$_{OH}$-O$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 3.5),
    )
    str_pair = 'cc_o_oh.h_w'
    run(
        dict_data = {
            'AIMD': f'{dir_aimd}rdf.{str_pair}.ave.csv',
            'DPMD': f'{dir_dpmd}rdf.{str_pair}.ave.csv',
        },
        str_text = r'O$_{OH}$-H$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 2.5),
    )
    str_pair = 'cc_o_c.h_w'
    run(
        dict_data = {
            'AIMD': f'{dir_aimd}rdf.{str_pair}.ave.csv',
            'DPMD': f'{dir_dpmd}rdf.{str_pair}.ave.csv',
        },
        str_text = r'$^=$O-H$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 2.5),
    )
    str_pair = 'o_w.o_w'
    run(
        dict_data = {
            'AIMD': f'{dir_aimd}rdf.{str_pair}.ave.csv',
            'DPMD': f'{dir_dpmd}rdf.{str_pair}.ave.csv',
        },
        str_text = r'O$_W$-O$_W$',
        str_pair = str_pair,
        tup_ylim = (0, 4),
    )

    plt.show()

main()

