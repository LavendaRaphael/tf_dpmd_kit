import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

with open('COLVAR', 'r') as colvar:
    list_header = colvar.readline().split()[2:]
data = np.genfromtxt("COLVAR", dtype=None, names=list_header)
print(data.dtype)
#data['time'] /= 2.0

list_label = [
    #"dist_pro_c",
    #'cn_o_126_h',
    #'cn_o_127_h',
    #'cn_o_128_h',
    #'del_cn_o_h',
    #'cn_c_pro',
    #'cn_c_pro_1',
    #'cost_o_126_h',
    #'cost_o_127_pro',
    #'cost_o_128_pro',
    #'cost_o_127_vh',
    #'cost_o_128_vh',
    #'cost_o_127_newp',
    #'cost_o_128_newp',
    #'restraintbias',
    #'uwallbias',
    #'uwallforce2',
    'dist_pro_o_127',
    #'dist_pro_o_128',
]

for str_label in list_label:
    ax.plot(data['time'], data[str_label], label=str_label)
#str_save = 'cv_t_com.pdf'

ax.legend()
ax.set_xlabel('Time(ps)')
ax.set_ylabel('CV')
#plt.savefig(str_save, bbox_inches='tight')
plt.show()
