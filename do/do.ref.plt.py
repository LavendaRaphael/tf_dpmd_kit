import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def def_pka(
    float_T: float
) -> float:
    return -373.604+16500.0/float_T + 56.478*np.log(float_T)

def def_deltag(
    float_T: float,
) -> float:
    float_Avogadro = 6.02214076e23
    # J*K^-1
    float_Kb = 1.380649e-23
    # KJ*K^-1*mol^-1
    float_R = float_Kb*float_Avogadro/1000.0
    return float_R * np.log(10) * (16500 - 373.604*float_T + 56.478*float_T*np.log(float_T))

np_T = np.linspace(
    start = 273+5,
    stop = 273+40,
    num = 100
)

matplotlib.rcParams['font.size']=15
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=["Arial"]

fig, ax = plt.subplots()

#ax.plot( np_T, def_pka(np_T), label='-373.604 + 16500/T + 56.478lnT')
#ax.set_ylabel('pKa')

ax.plot( np_T, def_deltag(np_T), label='kT ln10(16500 - 373.604T + 56.478TlnT)')
ax.set_ylabel(r'$\Delta$G (KJ/mol)')

ax.legend()
ax.set_xlabel('T (K)')
fig.savefig('ref.deltag.pdf', bbox_inches='tight')
plt.show()
