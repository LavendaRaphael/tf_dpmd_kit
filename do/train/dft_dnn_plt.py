from matplotlib import pyplot as plt
import tf_dpmd_kit.train

#'''
tf_dpmd_kit.train.dptest_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    str_save = 'dptest.e.pdf',
    tup_inset_xlim = (-1.1,1.1)
)
#'''
#'''
tf_dpmd_kit.train.dptest_plt(
    str_file = 'dptest.f.out',
    str_save = 'dptest.f.pdf',
)
#'''

plt.show()
