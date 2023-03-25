import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
import pandas as pd

plot.set_rcparam()
cm = 1/2.54

def run():

    df_data = pd.read_csv('carbonic_life_time.csv')

    fig, ax = plt.subplots()

    list_header = [ 'CT','TT', 'CC']
    for header in list_header:
        #ax.hist(df_data[header], bins='auto', density=True, label=header)
        ax.hist(df_data[header], bins=100, label=header, alpha=0.5, range=(0,9000))

    ax.legend()

    plot.save(
        fig,
        tup_size = (8.6*cm, 5*cm)
    )

run()

plt.show()

