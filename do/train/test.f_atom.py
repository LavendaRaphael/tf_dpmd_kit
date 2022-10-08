import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def def_plt(
) -> None:
    
    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    np_data = np.loadtxt('dptest.f.out')
    del_data_xyz = np_data[:,3:6] - np_data[:,0:3]
    del_data = np.linalg.norm( del_data_xyz, axis=1 )
    del_data *= 1000

    ax.plot(del_data)

    ax.set_ylabel(r'|F$_{DP}$-F$_{DFT}$| (meV/Å)')
    ax.set_xlabel('Atom index')

    # filter large
    for int_i in range(len(del_data)):
        if del_data[int_i] > 300:
            print(int_i, del_data[int_i])

def_plt()

plt.show()
