from matplotlib import pyplot as plt
import numpy as np
import os

# data
list2d_dirdata = []
#'''
list2d_dirdata.append(['rdf.c.o_w.0010000_0360000.npy', 'C-O_w' ])
list2d_dirdata.append(['rdf.o_1.h_w.0010000_0360000.npy', 'O_=-H_w' ])
list2d_dirdata.append(['rdf.o_0_2.h_w.0010000_0360000.npy', 'O_h-O_w' ])
list2d_dirdata.append(['rdf.h_0_1.o_w.0010000_0360000.npy', 'H_c-O_w' ])
list2d_dirdata.append(['rdf.o_w.o_w.0010000_0360000.npy', 'O_w-O_w' ])
#'''

# common

fig, ax = plt.subplots()

for list_dirdata in list2d_dirdata:
    array_rdf = np.load( os.path.join(list_dirdata[0]) )
    ax.plot(
        array_rdf[0], 
        array_rdf[1], 
        label = list_dirdata[1],
        )

ax.legend()
ax.set_xlabel('r (Å)')
ax.set_ylabel('RDF')
ax.set_xlim(None)
ax.set_ylim(None)
#fig.set_size_inches(8, 4)
fig.savefig('rdf.pdf', bbox_inches='tight')
plt.show()

