import numpy as np
import math

def def_T2KbT(
    float_T: float,
) -> float:

    float_Avogadro = 6.02214076e23
    # J*K^-1
    float_Kb = 1.380649e-23
    # KJ*mol^-1
    float_KbT = float_Kb*float_Avogadro*float_T/1000.0
    return float_KbT

def def_integral(
    str_file: str,
    float_T: float,
    float_xlow: float,
    float_xup: float,
) -> float:

    if float_xup == '1M':
        float_Avogadro = 6.02214076e23
        float_xup = 1e27/float_Avogadro
        float_xup = (float_xup/(4/3*math.pi))**(1.0/3.0)
        print(f'float_xup = {float_xup}')

    np_data = np.loadtxt(str_file)    
    np_data[:,1] -= min(np_data[:,1])
 
    float_KbT = def_T2KbT(float_T)
    
    np_data[:,1] = np.exp(-np_data[:,1]/float_KbT)
 
    np_indices = np.searchsorted(np_data[:,0], [float_xlow, float_xup])
    print('xrange = ', np_data[np_indices[0],0], np_data[np_indices[1],0])
    slice_range = slice(np_indices[0], np_indices[1])
    float_integral = np.trapz(
        np_data[slice_range, 1],
        x = np_data[slice_range, 0]
    )

    float_integral = -np.log(float_integral) * float_KbT
    return float_integral

def def_pka(
    float_deltag: float,
    float_T: float,
) -> None:

    float_KbT = def_T2KbT(float_T)
    return float_deltag/(float_KbT*math.log(10))

float_T = 350

float_g1 = def_integral(
    str_file = 'fes.dist_vp_c.grid',
    float_T = float_T,
    float_xlow = 1,
    float_xup = 2,
)
print('float_g1', float_g1)
float_g2 = def_integral(
    str_file = 'fes.dist_vp_c.grid',
    float_T = float_T,
    float_xlow = 2.2,
    float_xup = '1M',
)
print('float_g2', float_g2)
float_deltag = float_g2-float_g1
float_pka = def_pka(
    float_deltag = float_deltag,
    float_T = float_T,
)
print('float_deltag', float_deltag)
print('float_pka', float_pka)
