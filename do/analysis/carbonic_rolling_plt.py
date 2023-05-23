from tf_dpmd_kit import plot
from tf_dpmd_kit import analysis
import matplotlib.pyplot as plt

def run(ax, xscale, xlabel):

    analysis.carbonic_rolling_plt(
        ax,
        float_xscale = xscale,
        file_data = 'carbonic_state.csv',
        int_window = 100,
        list_header = ['CO3','HCO3','CC','CT','TT','H3CO3'],
        list_ypos = [0, 1, 2, 3, 4, 5],
    )
    ax.set_xlabel(xlabel)

def main():

    plot.set_rcparam()
    cm = 1/2.54
    
    fig, ax = plt.subplots(figsize=(8.6*cm, 3.5*cm))
    
    #run(ax, 0.0004837769, 'Time (ps)')
    #run(ax, 0.00001, 'Time (ns)')
    run(ax, 0.000005, 'Time (ns)')
    
    ax.set_yticklabels([r'CO$_3^{2-}$', 'HCO$_3^-$', 'CC','CT','TT', r'H$_3$CO$_3$'])
    
    plot.save(
        fig,
        file_save = 'carbonic_rolling',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()
