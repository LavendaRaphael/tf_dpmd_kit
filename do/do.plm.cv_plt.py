import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
'''
data = np.genfromtxt("COLVAR", dtype=None, names=['time','cv1','cv2','cv3'], usecols=(0,1,2,3))
data['time'] /= 2.0
ax.plot(data['time'], data['cv1'], label='cn(O_126-H)')
ax.plot(data['time'], data['cv2'], label='cn(O_127-H)')
ax.plot(data['time'], data['cv3'], label='cn(O_128-H)')
str_save = 'cv_cn.pdf'
#'''
'''
data = np.genfromtxt("COLVAR", dtype=None, names=['time','cv1','cv2','cv3'], usecols=(0,4,5,6))
data['time'] /= 2.0
ax.plot(data['time'], data['cv1'], label='T(O_126-H)')
ax.plot(data['time'], data['cv2'], label='T(O_127-H)')
ax.plot(data['time'], data['cv3'], label='T(O_128-H)')
str_save = 'cv_t.pdf'
#'''
#'''
data = np.genfromtxt("COLVAR", dtype=None, names=['time','cv1','cv2','cv3'], usecols=(0,7,8,9))
data['time'] /= 2.0
ax.plot(data['time'], data['cv1'], label='T(O_126-V)')
ax.plot(data['time'], data['cv2'], label='T(O_127-V)')
ax.plot(data['time'], data['cv3'], label='T(O_128-V)')
str_save = 'cv_t_com.pdf'
#'''

ax.legend()
ax.set_xlabel('Time(ps)')
ax.set_ylabel('CV')
plt.savefig(str_save, bbox_inches='tight')
plt.show()
