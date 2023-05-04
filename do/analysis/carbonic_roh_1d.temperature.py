from tf_dpmd_kit import plot
import matplotlib.pyplot as plt
import matplotlib as mpl

def run(ax):
    dict_color = plot.dict_color_temperature(
        tup_colormap = (290, 350),
        dict_temperature = {
            '290K': 290,
            '310K': 310,
            '330K': 330,
            '350K': 350,
        }
    )
    
    str_tmp = 'carbonic/carbonic_roh_1d.csv'
    plot.plt_compare(
        ax,
        dict_data = {
            '290K': f'../290K/{str_tmp}',
            '310K': f'../310K/{str_tmp}',
            '330K': f'../330K/{str_tmp}',
            '350K': f'../350K/{str_tmp}',
        },
        dict_color = dict_color,
    )

    ax.set_xlabel(r'R(OH) (Å)')
    ax.set_ylabel('Free energy (kJ/mol)')
    ax.legend()
    
def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        #file_save = 'carbonic_roh_1d',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
