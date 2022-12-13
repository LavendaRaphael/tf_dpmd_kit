import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def dptest_plt(
    str_file: str,
    int_natoms: int = None,
    str_save: str = None,
    str_plt_type: str = 'inset',
    tup_inset_xlim: tuple = None,
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

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

    fig, ax = plt.subplots()

    float_rmse_e, float_rmse_f, float_rmse_v = get_rmse()

    if (str_mode=='e'):
        # per atom
        np_data_new = np_data / int_natoms
        # eV to meV
        float_rmse_e *= 1000
        np_data_new *= 1000
        np_data_new -= np.average( np_data_new[:,0] )

        tup_dft_dnn_label = ('DFT energy (meV/atom)', 'DNN energy (meV/atom)')
        str_label = f'Energy RMSE = {float_rmse_e:.3f} meV/atom'
        str_inset_xlabel = 'error (meV/atom)'
        str_hist_xlabel = r'E$_{DP}$-E$_{DFT}$ (meV/atom)'

    elif (str_mode=='f'):
        np_data_new = np_data.reshape((np_data.shape[0]*3, 2), order='F')
        
        tup_dft_dnn_label = ('DFT force (eV/Å)', 'DNN force (eV/Å)')
        str_label = f'Force RMSE = {float_rmse_f:.3f} eV/Å'
        str_inset_xlabel = 'error (eV/Å)'
        str_hist_xlabel = r'|F$_{DP}$-F$_{DFT}$| (eV/Å)'

    del_data = np_data_new[:,1] - np_data_new[:,0]
    if (str_plt_type == 'dft_dnn' or str_plt_type == 'inset'):
        ax.plot(
            np_data_new[:,0],
            np_data_new[:,1],
            label =str_label,
            marker='o',
            linestyle = '',
            markersize=2,
        )
        ax.set_xlabel(tup_dft_dnn_label[0])
        ax.set_ylabel(tup_dft_dnn_label[1])
        ax.axline([0, 0], [1, 1], color='black', linestyle='--')
        if str_plt_type == 'inset':
            axin = ax.inset_axes((0.65,0.2,0.3,0.3))
            axin.hist(
                del_data,
                bins = 'auto',
                density = True
            )
            axin.set_xlim(tup_inset_xlim)
            axin.set_xlabel(str_inset_xlabel)
            axin.set_ylabel('Probability')

    elif (str_plt_type == 'hist'):
        ax.hist(
            del_data,
            label = str_label,
            bins = 'auto',
            density = True
        )
        ax.set_xlabel(str_hist_xlabel)
        ax.set_ylabel('Probability')
    else:
        raise

    ax.legend()
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')

def get_rmse():
    with open('log', 'r') as fp:
        for str_line in fp:
            if 'Energy RMSE/Natoms' in str_line:
                float_rmse_e = float(str_line.split()[-2])
            elif 'Force  RMSE' in str_line:
                float_rmse_f = float(str_line.split()[-2])
            elif 'Virial RMSE/Natoms' in str_line:
                float_rmse_v = float(str_line.split()[-2])
    return float_rmse_e, float_rmse_f, float_rmse_v

def def_plt_lcurve(
    str_save: str=None
) -> None:

    matplotlib.rcParams['font.size']=15
    matplotlib.rcParams['font.family']='sans-serif'
    matplotlib.rcParams['font.sans-serif']=["Arial"]

    fig, ax = plt.subplots()

    data = np.genfromtxt("lcurve.out", names=True)
    for name in ['rmse_e_trn','rmse_f_trn']:
        ax.scatter(data['step'], data[name], label=name, s=0.5)
    ax.legend(loc='lower left')
    ax.set_xlabel('Step')
    ax.set_ylabel('Loss')
    ax.set_yscale('log')
    ax.set_xlim((None,None))
    ax.grid()
    if str_save:
        fig.savefig(str_save, bbox_inches='tight')
