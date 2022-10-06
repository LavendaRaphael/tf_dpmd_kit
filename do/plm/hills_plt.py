import numpy as np
import matplotlib.pyplot as plt

fig, axe = plt.subplots()
with open('HILLS', 'r') as colvar:
    list_header = colvar.readline().split()[2:]
data = np.genfromtxt("HILLS", dtype=None, names=list_header, invalid_raise=False)
#print(data.dtype)
axe.scatter(data['time'], data["height"], s=0.5)
axe.set_xlabel('Time (ps)')
axe.set_ylabel('Height (kJ)')
#fig.savefig('plm.hills.pdf', bbox_inches='tight')

plt.show()
