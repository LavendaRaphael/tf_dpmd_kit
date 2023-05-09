from tf_dpmd_kit import plot
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

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
    
    for dir_x in ['290K', '310K', '330K', '350K']:
        df = pd.read_csv(f'../{dir_x}/carbonic/carbonic_roh_1d.csv')
        ax.plot(df['roh0(ang)'], df['H2CO3'], color=dict_color[dir_x], lw=1, label=dir_x)

    ax.set_xlabel(r'R(OH) (Å)')
    ax.set_ylabel('Free energy (kcal/mol)')
    ax.legend()
    
def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300

    fig, ax = plt.subplots(figsize = (8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_roh_1d.temperature',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
