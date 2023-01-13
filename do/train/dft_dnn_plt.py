from matplotlib import pyplot as plt
import tf_dpmd_kit.train

#'''
tf_dpmd_kit.train.dptest_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    str_save = 'dptest.e.dft_dnn.svg',
    str_plt_type = 'dft_dnn',
)
tf_dpmd_kit.train.dptest_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    str_save = 'dptest.e.hist.svg',
    str_plt_type = 'hist',
    tup_hist_xlim = (-1.1, 1.1)
)
#'''
#'''
tf_dpmd_kit.train.dptest_plt(
    str_file = 'dptest.f.out',
    str_save = 'dptest.f.dft_dnn.png',
    str_plt_type = 'dft_dnn',
)
tf_dpmd_kit.train.dptest_plt(
    str_file = 'dptest.f.out',
    str_save = 'dptest.f.hist.svg',
    str_plt_type = 'hist',
)
#'''

plt.show()
