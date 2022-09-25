import numpy as np

def def_integral(
    str_file: str,
    float_T: float,
    float_xlow: float,
    float_xup: float,
) -> None:
    np_data = np.loadtxt(str_file)
    
    np_data[:,1] -= min(np_data[:,1])
    
    float_Avogadro = 6.02214076e23
    # J*K^-1
    float_Kb = 1.380649e-23
    # KJ*mol^-1
    float_KbT = float_Kb*float_Avogadro*float_T/1000.0
    
    np_data[:,1] = np.exp(-np_data[:,1]/float_KbT)
 
    np_indices = np.searchsorted(np_data[:,0], [float_xlow, float_xup])
    slice_range = slice(np_indices[0], np_indices[1])
    float_integral = np.trapz(
        np_data[slice_range, 1],
        x = np_data[slice_range, 0]
    )

    float_integral = -np.log(float_integral) * float_KbT
    print(float_integral)
    return float_integral

def_integral(
    str_file = 'fes.dist_vp_c.grid',
    float_T = 370,
    float_xlow = 1,
    float_xup = 2,
)
def_integral(
    str_file = 'fes.dist_vp_c.grid',
    float_T = 370,
    float_xlow = 2.2,
    float_xup = 7,
)

