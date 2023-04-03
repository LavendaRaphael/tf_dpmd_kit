from matplotlib import pyplot as plt
from tf_dpmd_kit import train
from tf_dpmd_kit import plot

plot.set_rcparam()
cm = 1/2.54

#'''
fig, ax = plt.subplots(figsize=(4.3*cm, 4*cm))
train.dptest_parity_plt(
    ax,
    str_file = 'dptest.e.out',
    int_natoms = 384,
    float_lw = 0.75,
    list_ticks = [-5, 0, 5],
)
plot.save(
    fig,
    file_save = 'dptest.e.parity',
    list_type = ['pdf', 'svg']
)

fig, ax = plt.subplots(figsize=(4.3*cm, 4*cm))
train.dptest_parity_plt(
    ax,
    str_file = 'dptest.f.out',
    float_lw = 0.75,
    list_ticks = [-5, 0, 5],
)
plot.save(
    fig,
    file_save = 'dptest.f.parity',
    list_type = ['pdf', 'svg']
)
#'''

plt.show()
