import numpy as np
import math
from tf_dpmd_kit import plm

def fes_integral(
    str_file: str,
    float_xlow: float,
    float_xup: float,
    float_T: float = None,
    str_plmin: str = '../../plm.in',
    str_plmlog: str = '../../plm.log'
) -> float:

    if float_xup == '1M':
        float_Avogadro = 6.02214076e23
        float_xup = 1e27/float_Avogadro
        float_xup = (float_xup/(4/3*math.pi))**(1.0/3.0)
        print(f'float_xup = {float_xup}')

    np_data = np.loadtxt(str_file)    
    np_data[:,1] -= min(np_data[:,1])

    if not float_T:
        float_T, float_KbT = plm.get_temperature( str_plmin, str_plmlog )
    else:
        float_KbT = plm.T2KbT(float_T)
    
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

float_g1 = fes_integral(
    str_file = 'fes.dist_vp_c.grid',
    float_xlow = 1.16,
    float_xup = 1.66,
)
print('float_g1', float_g1)

float_g2 = fes_integral(
    str_file = 'fes.dist_vp_c.grid',
    float_xlow = 2.56,
    float_xup = '1M',
)
print('float_g2', float_g2)

float_deltag = float_g2-float_g1
float_pka = plm.deltag_to_pka(
    float_deltag = float_deltag,
)
print('float_deltag', float_deltag)
print('float_pka', float_pka)
