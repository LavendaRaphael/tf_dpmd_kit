from tf_dpmd_kit import analysis
from tf_dpmd_kit import plot
import matplotlib.pyplot as plt

def sub(
    ax,
):

    df = analysis.carbonic_survival(
        list_file = [
            '../CC/carbonic/carbonic_lifedata.csv',
            '../CT/carbonic/carbonic_lifedata.csv',
            '../TT/carbonic/carbonic_lifedata.csv',
        ]
    )

    analysis.carbonic_lifetime(df)

    list_state = ['CC', 'CT', 'TT','H2CO3', 'HCO3']
    dict_color = {
        'CC': 'tab:blue',
        'CT': 'tab:orange',
        'TT': 'tab:green',
        'H2CO3': 'tab:red',
        'HCO3': 'tab:purple',
    }
    dict_label = {
        'H2CO3': r'H$_2$CO$_3$',
        'HCO3': r'HCO$_3^-$',
    }
    analysis.carbonic_survival_plt(
        ax = ax,
        df = df,
        list_state = list_state,
        dict_color = dict_color,
        dict_label = dict_label,
    )

def run():

    plot.set_rcparam()
    cm = 1/2.54

    fig, ax = plt.subplots(figsize=(8.6*cm, 5*cm))

    sub(ax)

    plot.save(
        fig,
        file_save = 'carbonic_survival',
        list_type = ['pdf', 'svg']
    )
    plt.show()

run()
