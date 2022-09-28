import numpy as np
from matplotlib import pyplot as plt

def def_f(nn, mm, r0, r):
    return (1 - np.power(r/r0, nn)) / (1 - np.power(r/r0, mm))

r = np.linspace(
    start = 0.5,
    stop = 2,
    num = 100
)

fig, ax = plt.subplots()

ax.plot( r, def_f(6, 12, 1.3, r), label='NN, MM = 6, 12', linewidth=2)
ax.plot( r, def_f(12, 24, 1.3, r), label='NN, MM = 12, 24', linewidth=2)

ax.set_xlabel('R')
ax.set_ylabel('s(R)')
ax.legend()
fig.savefig('plm.para.cn.pdf', bbox_inches='tight')
plt.show()

