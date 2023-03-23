import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
import pandas as pd

plot.set_rcparam()
cm = 1/2.54

def run():

    pd_data = pd.read_csv('carbonic_lifetime.csv')

    fig, ax = plt.subplots()

    ax.hist(pd_data['TT'], bins='auto')
    ax.hist(pd_data['CT'], bins='auto')
    ax.hist(pd_data['CC'], bins='auto')

    plot.save(
        fig,
        tup_size = (8.6*cm, 5*cm)
    )

run()

plt.show()

