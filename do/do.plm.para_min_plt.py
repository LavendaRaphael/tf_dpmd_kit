import numpy as np
from matplotlib import pyplot as plt

def def_f(n1, n2, beta, cn):
    return beta/(np.log( n1*np.exp(beta/cn) + n2*np.exp(beta/2) ))

def def_g(n1, n2, beta, cn):
    return -1.0/beta*np.log( n1*np.exp(-beta*cn) + n2*np.exp(-beta*2) )

cn = np.linspace(
    start = 1,
    stop = 3,
    num = 100
)
n1 = 1
n2 = 125

fig, ax = plt.subplots()

ax.plot( cn, def_f(n1, n2, 50, cn), label='beta = 50', linewidth=2)
ax.plot( cn, def_f(n1, n2, 100, cn), label='beta = 100', linewidth=2)
ax.plot( cn, def_f(n1, n2, 200, cn), label='beta = 200', linewidth=2)

#ax.plot( cn, def_g(n1, n2, 10, cn), label='g 10', linewidth=2)    
#ax.plot( cn, def_g(n1, n2, 50, cn), label='g 50', linewidth=2)  
#ax.plot( cn, def_g(n1, n2, 100, cn), label='g 100', linewidth=2)

ax.set_xlabel('s0')
ax.set_ylabel('min')
ax.legend()
fig.savefig('plm.para.min.pdf', bbox_inches='tight')
plt.show()

