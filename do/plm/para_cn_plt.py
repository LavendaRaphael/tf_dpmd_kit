import numpy as np
from matplotlib import pyplot as plt

def def_f(nn, mm, r0, r):
    return (1 - np.power(r/r0, nn)) / (1 - np.power(r/r0, mm))

def cn_plt(
    list_mn: list
) -> None:
    r = np.linspace(
        start = 0.5,
        stop = 2,
        num = 100
    )
    
    fig, ax = plt.subplots()
    
    for tup_mn in list_mn:
        mm = tup_mn[0]
        nn = tup_mn[1]
        r0 = tup_mn[2]
        ax.plot( r, def_f(mm, nn, r0, r), label=f'NN, MM, r0 = {mm}, {nn}, {r0}', linewidth=2)
    
    ax.set_xlabel('R')
    ax.set_ylabel('s(R)')
    ax.legend()
    fig.savefig('plm.para.cn.pdf', bbox_inches='tight')

cn_plt(
    list_mn = [
        (6,18,1.4),
        (6,12,1.4),
    ]
)
plt.show()


