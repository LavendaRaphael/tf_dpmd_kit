from tf_dpmd_kit import train
from tf_dpmd_kit import plot
from matplotlib import pyplot as plt

plot.set_rcparam()

fig, ax = train.plt_lcurve()

cm = 1/2.54
plot.save(
    fig,
    tup_size = (8.6*cm, 5*cm),
    str_save = 'lcurve.png',
)

plt.show()
