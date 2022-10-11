import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def check_locate(
    int_id: int,
):
    list_dpnpy = [
        ('      data.init/H2CO3_CC_H2O_126/dpnpy', 237), 
        ('      data.init/H2CO3_CT_H2O_126/dpnpy', 354), 
        ('      data.init/H2CO3_TT_H2O_126/dpnpy',  61), 
        ('    data.init/H2O_128/dpnpy/iter.0.004',  13), 
        ('    data.init/H2O_128/dpnpy/iter.1.005', 300), 
        ('    data.init/H2O_128/dpnpy/iter.2.006', 236), 
        ('    data.init/H2O_128/dpnpy/iter.3.007', 107), 
        ('    data.init/H2O_128/dpnpy/iter.4.009', 149), 
        ('    data.init/H2O_128/dpnpy/iter.5.011',  20), 
        ('    data.init/H2O_128/dpnpy/iter.6.013',   1), 
        ('    data.init/H3O_OH_far_H2O_126/dpnpy', 366), 
        ('   data.init/H3O_OH_near_H2O_126/dpnpy', 188), 
        ('  data.init/HCO3_H3O_far_H2O_125/dpnpy',  99), 
        (' data.init/HCO3_H3O_near_H2O_125/dpnpy', 386), 
        (' data.iters/iter.000000/02.fp/data.000', 300), 
        (' data.iters/iter.000000/02.fp/data.004', 300), 
        (' data.iters/iter.000000/02.fp/data.008', 300), 
        (' data.iters/iter.000000/02.fp/data.012',  15), 
        (' data.iters/iter.000000/02.fp/data.016', 300), 
        (' data.iters/iter.000000/02.fp/data.020', 300), 
        (' data.iters/iter.000001/02.fp/data.001', 300), 
        (' data.iters/iter.000001/02.fp/data.005', 300), 
        (' data.iters/iter.000001/02.fp/data.009', 300), 
        (' data.iters/iter.000001/02.fp/data.013', 225), 
        (' data.iters/iter.000001/02.fp/data.017', 300), 
        (' data.iters/iter.000001/02.fp/data.021', 300), 
        (' data.iters/iter.000002/02.fp/data.002', 300), 
        (' data.iters/iter.000002/02.fp/data.006', 300), 
        (' data.iters/iter.000002/02.fp/data.010', 300), 
        (' data.iters/iter.000002/02.fp/data.014',  27), 
        (' data.iters/iter.000002/02.fp/data.018', 300), 
        (' data.iters/iter.000002/02.fp/data.022', 300), 
        (' data.iters/iter.000003/02.fp/data.003', 300), 
        (' data.iters/iter.000003/02.fp/data.007', 300), 
        (' data.iters/iter.000003/02.fp/data.011', 300), 
        (' data.iters/iter.000003/02.fp/data.015', 300), 
        (' data.iters/iter.000003/02.fp/data.019', 300), 
        (' data.iters/iter.000003/02.fp/data.023', 300), 
        (' data.iters/iter.000004/02.fp/data.024', 300), 
        (' data.iters/iter.000004/02.fp/data.028', 300), 
        (' data.iters/iter.000004/02.fp/data.032', 300), 
        (' data.iters/iter.000004/02.fp/data.036',  79), 
        (' data.iters/iter.000004/02.fp/data.040', 300), 
        (' data.iters/iter.000004/02.fp/data.044', 300), 
        (' data.iters/iter.000005/02.fp/data.025', 300), 
        (' data.iters/iter.000005/02.fp/data.029', 300), 
        (' data.iters/iter.000005/02.fp/data.033', 300), 
        (' data.iters/iter.000005/02.fp/data.037',  52), 
        (' data.iters/iter.000005/02.fp/data.041', 300), 
        (' data.iters/iter.000005/02.fp/data.045', 300), 
        (' data.iters/iter.000006/02.fp/data.026', 222), 
        (' data.iters/iter.000006/02.fp/data.030', 148), 
        (' data.iters/iter.000006/02.fp/data.034', 141), 
        (' data.iters/iter.000006/02.fp/data.038',  46), 
        (' data.iters/iter.000006/02.fp/data.042', 230), 
        (' data.iters/iter.000006/02.fp/data.046', 206), 
        (' data.iters/iter.000007/02.fp/data.027', 300), 
        (' data.iters/iter.000007/02.fp/data.031', 300), 
        (' data.iters/iter.000007/02.fp/data.035', 300), 
        (' data.iters/iter.000007/02.fp/data.039', 167), 
        (' data.iters/iter.000007/02.fp/data.043', 300), 
        (' data.iters/iter.000007/02.fp/data.047', 300), 
        (' data.iters/iter.000008/02.fp/data.048', 279), 
        (' data.iters/iter.000008/02.fp/data.052', 279), 
        (' data.iters/iter.000008/02.fp/data.056', 183), 
        (' data.iters/iter.000009/02.fp/data.048', 300), 
        (' data.iters/iter.000009/02.fp/data.052', 158), 
        (' data.iters/iter.000009/02.fp/data.056', 124), 
        (' data.iters/iter.000009/02.fp/data.060', 109), 
        (' data.iters/iter.000010/02.fp/data.049', 115), 
        (' data.iters/iter.000010/02.fp/data.053',  86), 
        (' data.iters/iter.000010/02.fp/data.057',  65), 
        (' data.iters/iter.000010/02.fp/data.061', 103), 
        (' data.iters/iter.000011/02.fp/data.050', 243), 
        (' data.iters/iter.000011/02.fp/data.054', 215), 
        (' data.iters/iter.000011/02.fp/data.058', 210), 
        (' data.iters/iter.000011/02.fp/data.062', 172), 
        (' data.iters/iter.000012/02.fp/data.051', 290), 
        (' data.iters/iter.000012/02.fp/data.055', 253), 
        (' data.iters/iter.000012/02.fp/data.059', 245), 
        (' data.iters/iter.000012/02.fp/data.063', 206), 
        (' data.iters/iter.000013/02.fp/data.064',  84), 
        (' data.iters/iter.000013/02.fp/data.068', 150), 
        (' data.iters/iter.000013/02.fp/data.072',  68), 
        (' data.iters/iter.000013/02.fp/data.076',  91), 
        (' data.iters/iter.000014/02.fp/data.065', 300), 
        (' data.iters/iter.000014/02.fp/data.069', 298), 
        (' data.iters/iter.000014/02.fp/data.073', 298), 
        (' data.iters/iter.000014/02.fp/data.077', 300), 
        (' data.iters/iter.000015/02.fp/data.066', 298), 
        (' data.iters/iter.000015/02.fp/data.070', 300), 
        (' data.iters/iter.000015/02.fp/data.074', 299), 
        (' data.iters/iter.000015/02.fp/data.078', 299), 
        (' data.iters/iter.000016/02.fp/data.067', 300), 
        (' data.iters/iter.000016/02.fp/data.071', 300), 
        (' data.iters/iter.000016/02.fp/data.075', 300), 
        (' data.iters/iter.000016/02.fp/data.079', 300), 
        (' data.iters/iter.000017/02.fp/data.067',  10), 
        (' data.iters/iter.000017/02.fp/data.071',  13), 
        (' data.iters/iter.000017/02.fp/data.075',  10), 
        (' data.iters/iter.000017/02.fp/data.079',  20), 
    ]
    int_sum = 0
    for tup_dpnpy in list_dpnpy:
        int_sum += tup_dpnpy[1]

    for tup_dpnpy in list_dpnpy:
        if int_id < tup_dpnpy[1]:
            return tup_dpnpy[0], int_id
        else:
            int_id -= tup_dpnpy[1]
    raise

def def_plt(
    float_rmse: float,
    str_file: str,
    int_natoms: int,
    str_save: str = None,
    tup_xlim: tuple = None,
    tup_ylim: tuple = None,
    int_bins: int = 'auto',
) -> None:
    
    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    with open(str_file, 'r') as file_open:
        list_line = file_open.readline().split()
        if list_line[-1] == 'pred_e':
            str_mode = 'e'
        elif list_line[-1] == 'pred_fz':
            str_mode = 'f'
        else:
            print(list_line)
            raise
    np_data = np.loadtxt(str_file)
    if (str_mode=='e'):
        # per atom
        np_data /= int_natoms
        # eV to meV
        np_data *= 1000
        del_data = np_data[:,1] - np_data[:,0]
        ax.set_xlabel( r'E$_{DP}$-E$_{DFT}$ (meV/atom)')
        str_label = f'Energy RMSE = {float_rmse:.3f} meV/atom'
        for int_i in range(len(del_data)):
            if abs(del_data[int_i]) > 1.5:
                str_locate, int_locate = check_locate(int_i)
                print(f'{str_locate} {int_locate} {del_data[int_i]:.1f}')
    elif (str_mode=='f'):
        del_data_xyz = np_data[:,3:6] - np_data[:,0:3]
        del_data = np.linalg.norm( del_data_xyz, axis=1 )
        del_data *= 1000
        ax.set_xlabel(r'|F$_{DP}$-F$_{DFT}$| (meV/Å)')
        str_label = f'Force RMSE = {float_rmse:.1f} meV/Å'
        for int_i in range(len(del_data)):
            if del_data[int_i] > 1000:
                int_id = int_i // int_natoms
                str_locate, int_locate = check_locate(int_id)
                print(f'{str_locate} {int_locate} {del_data[int_i]:.1f}')

    #float_std = np.std(del_data)
    #float_mean = np.mean(del_data) 
    ax.hist(
        del_data,
        label = str_label,
        bins = int_bins,
        #density = True
    )

    ax.legend()
    ax.set_ylabel('Probabilty Density')
    ax.set_xlim(tup_xlim)
    ax.set_ylim(tup_ylim)
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

#'''
def_plt(
    str_file = 'dptest.e.out',
    int_natoms = 384,
    float_rmse = 0.440,
    #tup_xlim = (-1.5,1.5),
    tup_ylim = (0,10),
    #str_save = 'dptest.e.pdf'
)
#'''
#'''
def_plt(
    str_file = 'dptest.f.out',
    float_rmse = 60.0,
    int_natoms = 384,
    #tup_xlim = (0,250),
    tup_ylim = (0,10),
    #str_save = 'dptest.f.pdf'
)
#'''

plt.show()
