import numpy as np
import matplotlib.pyplot as plt

def def_plt_lcurve():
    fig, ax = plt.subplots()

    data = np.genfromtxt("lcurve.out", names=True)
    for name in data.dtype.names[2:-2]:
        ax.plot(data['step'], data[name], label=name)
    ax.legend()
    ax.set_xlabel('Step')
    ax.set_ylabel('Loss')
    ax.set_xscale('symlog')
    ax.set_yscale('log')
    ax.set_xlim((1,None))
    ax.grid()
    fig.savefig('lcurve.pdf', bbox_inches='tight')
    plt.show()
