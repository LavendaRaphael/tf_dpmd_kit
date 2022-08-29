import numpy as np

np_data = np.genfromtxt("rdf.c.o_w.csv", names=True)

print(np_data.dtype)
