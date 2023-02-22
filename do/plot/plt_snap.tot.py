import matplotlib.pyplot as plt
from tf_dpmd_kit import plot
import matplotlib as mpl
import numpy as np

plot.set_rcparam()
mpl.rcParams['figure.dpi'] = 300
fig, ax = plt.subplots()
ax.axis('off')

dir_cp  = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_TT_H2O_126/plm/'
dir_tt  = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_TT_H2O_126/330K/plm/'
dir_tt0 = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_TT_H2O_126.0/330K/plm/'
dir_ct  = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_CT_H2O_126/330K/plm/'
dir_ct0 = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_CT_H2O_126.0/330K/plm/'
dir_cc  = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_CC_H2O_126/330K/plm/'
dir_cc0 = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md/H2CO3_CC_H2O_126.0/330K/plm/'

dy = 0.25
dx = 1/7

y0 = 1-dy/2
y1 = 0.5
y2 = dy/2
y3 = y0-0.03
y4 = y1-0.03

pos_00 = np.array([0.5*dx, y0])
pos_01 = np.array([1.5*dx, y3])
pos_02 = np.array([3  *dx, y0])
pos_03 = np.array([4.5*dx, y3])
pos_04 = np.array([5.5*dx, y0])
pos_05 = np.array([6.5*dx, y0])

pos_10 = np.array([0.5*dx, y1])
pos_11 = np.array([1.5*dx, y4])
pos_12 = np.array([2.5*dx, y1])
pos_13 = np.array([3.5*dx, y1])
pos_14 = np.array([4.5*dx, y4])
pos_15 = np.array([5.5*dx, y1])

pos_20 = np.array([1.5*dx, y2])
pos_21 = np.array([   0.5, y2])

plot.inset_img(
    ax,
    dict_img = {
        dir_cc +'1.100001.png': ( pos_00[0]-0.5*dx, pos_00[1]-0.5*dy, dx, dy),
        dir_ct +'0.003523.png': ( pos_01[0]-0.5*dx, pos_01[1]-0.5*dy, dx, dy),
        dir_cp +   '60281.png': ( pos_02[0]-0.5*dx, pos_02[1]-0.5*dy, dx, dy),
        dir_cc +'0.004450.png': ( pos_03[0]-0.5*dx, pos_03[1]-0.5*dy, dx, dy),
        dir_cc +'0.003922.png': ( pos_04[0]-0.5*dx, pos_04[1]-0.5*dy, dx, dy),
        dir_cc +'0.003576.png': ( pos_05[0]-0.5*dx, pos_05[1]-0.5*dy, dx, dy),
        dir_cp +   '26878.png': ( pos_10[0]-0.5*dx, pos_10[1]-0.5*dy, dx, dy),
        dir_cp +   '59864.png': ( pos_11[0]-0.5*dx, pos_11[1]-0.5*dy, dx, dy),
        dir_ct0+'0.016090.png': ( pos_12[0]-0.5*dx, pos_12[1]-0.5*dy, dx, dy),
        dir_cp +   '60184.png': ( pos_13[0]-0.5*dx, pos_13[1]-0.5*dy, dx, dy),
        dir_cc +'0.004441.png': ( pos_14[0]-0.5*dx, pos_14[1]-0.5*dy, dx, dy),
        dir_cc +'0.004420.png': ( pos_15[0]-0.5*dx, pos_15[1]-0.5*dy, dx, dy),
        dir_cp +   '60075.png': ( pos_20[0]-0.5*dx, pos_20[1]-0.5*dy, dx, dy),
        dir_ct0+'0.079593.png': ( pos_21[0]-0.5*dx, pos_21[1]-0.5*dy, dx, dy),
    },
    dict_imgcolor = {
        dir_cc +'1.100001.png': 'tab:blue',
        dir_cp +   '60281.png': 'tab:orange',
        dir_cc +'0.003922.png': 'tab:green',
        dir_cc +'0.003576.png': 'tab:green',
    },
    bool_rot90 = True,
)

cm = 1/2.54
float_size_x = 8.6*2*cm
float_size_y = 10*cm

dx1 = 1000/1600 * float_size_y/float_size_x * dy
l = np.array([-dx1/2,    0])
r = np.array([ dx1/2,    0])
t = np.array([     0, dy/2])
b = np.array([     0,-dy/2])

str_arrowstyle = '<|-|>'
plot.add_arrow(
    ax,
    dict_arrow = {
        '01_00': [pos_01+l, pos_00+r],
        '02_01': [pos_02+l, pos_01+r],
        '03_02': [pos_03+l, pos_02+r],
        '05_04': [pos_05+l, pos_04+r],

        '11_10': [pos_11+l, pos_10+r],
        '12_11': [pos_12+l, pos_11+r],
        '14_13': [pos_14+l, pos_13+r],
        '15_14': [pos_15+l, pos_14+r],

        '10_00': [pos_10+t, pos_00+b],
        '11_01': [pos_11+t, pos_01+b],
        '12_02': [pos_12+t, pos_02+b],
        '13_02': [pos_13+t, pos_02+b],
        '14_03': [pos_14+t, pos_03+b],
        '15_04': [pos_15+t, pos_04+b],

        '10_21': [pos_10+b, pos_21+l+[0,0.1]],
        '11_20': [pos_11+b, pos_20+t],
        '11_21': [pos_11+b, pos_21+l+[0,0.12]],
        '20_12': [pos_20+t, pos_12+b],
        '12_21': [pos_12+b, pos_21+t],
        '13_21': [pos_13+b, pos_21+t],
        '14_21': [pos_14+b, pos_21+t],
        '21_15': [pos_21+r+[0,0.1], pos_15+b],
    },
    dict_arrowstyle = {
        '05_04': str_arrowstyle,
        '11_10': str_arrowstyle,
        '12_11': str_arrowstyle,
        '14_13': str_arrowstyle,
        '15_14': str_arrowstyle,
        '10_00': str_arrowstyle,
        '12_02': str_arrowstyle,
        '13_02': str_arrowstyle,
        '15_04': str_arrowstyle,
        '10_21': str_arrowstyle,
        '11_21': str_arrowstyle,
        '12_21': str_arrowstyle,
        '13_21': str_arrowstyle,
        '14_21': str_arrowstyle,
    }
)

plot.save(
    fig,
    str_save = 'tot.svg',
    tup_size = (float_size_x, float_size_y),
)
plot.save(
    fig,
    str_save = 'tot.pdf',
    tup_size = (float_size_x, float_size_y),
)
plt.show()
