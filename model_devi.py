from ase import Atoms
from deepmd.calculator import DP
from numpy import linalg as LA
from matplotlib import pyplot as plt 
import numpy as np
import ase.io
import dpdata

def def_plt_devi_f():
    data = np.genfromtxt("model_devi.out", names=True)
    for name in ["max_devi_f"]:
        plt.plot(data['step'], data[name])
    plt.legend()
    plt.xlabel('Step')
    plt.ylabel('max_devi_f')
    plt.grid()
    #plt.savefig('max_devi_f.pdf', bbox_inches='tight')
    plt.show()

def def_model_devi_atom(
        list_model,
        dp_sys,
        float_lowthred,
        str_save,
        float_ylim = None
        ):
    
    ase_atoms = dp_sys.to('ase/structure')[0]
    
    int_natom = len(ase_atoms.get_atomic_numbers())
    
    int_nmodel = len(list_model)

    # calculate force
    array3d_forces = np.zeros (shape = (int_nmodel,int_natom,3))
    for int_i in range(int_nmodel):
        ase_atoms.calc = DP( model=list_model[int_i] )
        array3d_forces[int_i] = ase_atoms.get_forces()
    
    # average
    array2d_forces_ave = np.average (array3d_forces, axis=0)
    
    # devi
    array3d_devi_f = np.zeros(shape=(int_nmodel,int_natom,3))
    for int_i in range(int_nmodel):
        array3d_devi_f[int_i] = array3d_forces[int_i] - array2d_forces_ave
    
    # norm square
    array2d_devi_f_norm = LA.norm(array3d_devi_f, axis=2)**2
    
    # average, final
    array1d_devi_f_final = np.sqrt(np.average (array2d_devi_f_norm, axis=0))
    
    # filter large
    list_devi_f = []
    for int_i in range(int_natom):
        if array1d_devi_f_final[int_i] > float_lowthred:
            list_devi_f.append( [int_i, array1d_devi_f_final[int_i]] )
    print(list_devi_f)    

    # plot
    x = np.arange(int_natom)
    fig, ax = plt.subplots()
    ax.plot(x, array1d_devi_f_final)
    ax.set_xlabel("Atom index")
    ax.set_ylabel("Devi_f")
    ax.set_ylim(top=float_ylim)
    plt.savefig(str_save, bbox_inches='tight')
    plt.show()
