from matplotlib import pyplot as plt
from tf_dpmd_kit import train
from tf_dpmd_kit import plot

plot.set_rcparam()

fig = train.dptest_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    str_plt_type = 'dft_dnn',
    float_lw = 0.75,
)
cm = 1/2.54
plot.save(
    fig,
    tup_size = (4.3*cm, 3.7*cm),
    str_save = 'dptest.e.dft_dnn.pdf',
)

fig = train.dptest_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    str_plt_type = 'hist',
    tup_hist_xlim = (-1.1, 1.1),
    float_lw = 0.75,
)
plot.save(
    fig,
    tup_size = (4.3*cm, 3.7*cm),
    str_save = 'dptest.e.hist.pdf',
)

fig = train.dptest_plt(
    str_file = 'dptest.f.out',
    str_plt_type = 'dft_dnn',
    float_lw = 0.75,
)
plot.save(
    fig,
    tup_size = (4.3*cm, 3.7*cm),
    str_save = 'dptest.f.dft_dnn.png',
)

fig = train.dptest_plt(
    str_file = 'dptest.f.out',
    str_plt_type = 'hist',
)
plot.save(
    fig,
    tup_size = (4.3*cm, 3.7*cm),
    str_save = 'dptest.f.hist.pdf',
)

plt.show()
