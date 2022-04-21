import numpy as np
import matplotlib.pyplot as plt

def def_plt_lcurve():
    data = np.genfromtxt("lcurve.out", names=True)
    for name in data.dtype.names[1:-1]:
        plt.plot(data['step'], data[name], label=name)
    plt.legend()
    plt.xlabel('Step')
    plt.ylabel('Loss')
    plt.xscale('symlog')
    plt.yscale('symlog')
    plt.grid()
    plt.savefig('lcurve.pdf', bbox_inches='tight')
    plt.show()
