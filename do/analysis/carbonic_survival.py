from tf_dpmd_kit import analysis
from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

def run(
    ax,
):
    list_state = ['TT', 'CT', 'CC','H2CO3', 'HCO3']
    dict_color = {
        'TT': 'tab:blue',
        'CT': 'tab:orange',
        'CC': 'tab:green',
        'H2CO3': 'tab:red',
        'HCO3': 'tab:purple',
    }
    dict_label = {
        'H2CO3': r'H$_2$CO$_3$',
        'HCO3': r'HCO$_3^-$',
    }

    analysis.carbonic_survival(
        ax = ax,
        list_state = list_state,
        dict_color = dict_color,
        dict_label = dict_label,
    )
def main():

    plot.set_rcparam()
    cm = 1/2.54
    
    fig, ax = plt.subplots(figsize=(8.6*cm, 5*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_survival',
        list_type = ['pdf', 'svg']
    )
    plt.show()

main()

