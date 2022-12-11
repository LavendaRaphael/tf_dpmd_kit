import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
import matplotlib 

def plot(
    dict_data: dict,
    str_save: str = None,
    float_scale: float = 1.0,
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]
    
    fig, ax = plt.subplots()
    
    for str_label, str_data in dict_data.items():
        np_data = np.loadtxt(str_data)
        ax.scatter( np_data[:,0]*float_scale, np_data[:,1], s=0.5, label=str_label)
    
    ax.set_xlabel('t (ps)')
    ax.set_ylabel('n H-bonds')

    if str_save:
        #fig.set_size_inches(7, 6)
        fig.savefig(str_save, bbox_inches='tight')

plot(
    dict_data = {
        '': 'hbonds.h_0_1.o_w.0000000_0057877.csv',
    }
    #str_save = 'hbonds.pdf',
    float_scale = 1/0.005
)

plt.show()

