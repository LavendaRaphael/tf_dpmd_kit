import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt
import matplotlib 

def plot(
    str_data: str,
    str_save: str = None,
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]
    
    fig, ax = plt.subplots()
    
    np_data = np.loadtxt(str_data)

    ax.plot( np_data[:,0], np_data[:,1], label = 'MSD', linewidth=2)
    
    '''
    np_time = np.np([200,600])
    np_index = np.np( np_time / (np_data[0][1]-np_data[0][0]), dtype=int )
    np_x = np_data[0][ np_index[0]: np_index[1] ]
    np_y = np_data[1][ np_index[0]: np_index[1] ]
    linear_model = linregress( np_x, np_y )
    float_diffussioncoef = linear_model.slope / 6.0
    ax.plot( np_x, linear_model.intercept + linear_model.slope*np_x , '--', linewidth=2, label = f'D = {float_diffussioncoef:.3f} Å$^2$/ps')
    '''
    ax.legend()
    ax.set_xlabel('t (ps)')
    ax.set_ylabel('Mean Square Displacement (Å$^2$)')

    if str_save:
        #fig.set_size_inches(7, 6)
        fig.savefig(str_save, bbox_inches='tight')

plot(
    str_data = 'msd.csv',
    str_save = 'msd.pdf'
)

plt.show()

