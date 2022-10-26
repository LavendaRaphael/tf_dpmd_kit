import numpy as np
import math
from tf_dpmd_kit import plm
import glob

def fes_integral(
    str_file: str,
    tup_xrange: tuple,
    float_T: float = None,
    str_in: str = '../plm.in',
    str_log: str = '../plm.log'
) -> float:


    float_xlow = tup_xrange[0]
    float_xup = tup_xrange[1]
    if float_xup == '1M':
        float_Avogadro = 6.02214076e23
        float_xup = 1e27/float_Avogadro
        float_xup = (float_xup/(4/3*math.pi))**(1.0/3.0)

    np_data = np.loadtxt(str_file)    
    np_data[:,1] -= min(np_data[:,1])

    if not float_T:
        float_T, float_KbT = plm.get_temperature( str_in, str_log )
    else:
        float_KbT = plm.T2KbT(float_T)
    
    np_data[:,1] = np.exp(-np_data[:,1]/float_KbT)
 
    np_indices = np.searchsorted(np_data[:,0], [float_xlow, float_xup])
    slice_range = slice(np_indices[0], np_indices[1])
    float_integral = np.trapz(
        np_data[slice_range, 1],
        x = np_data[slice_range, 0]
    )

    float_integral = -np.log(float_integral) * float_KbT
    return float_integral

def get_pka(
    str_file: str,
    tup_xrange1: tuple,
    tup_xrange2: tuple,
    float_T: float = None,
    str_in: str = '../plm.in',
    str_log: str = '../plm.log'
) -> (float, float):

    if not float_T:
        float_T, float_KbT = plm.get_temperature( str_in, str_log )

    float_g1 = fes_integral(
        str_file = str_file,
        tup_xrange = tup_xrange1,
        float_T = float_T,
        str_in = str_in,
        str_log = str_log
    )
    
    float_g2 = fes_integral(
        str_file = str_file,
        tup_xrange = tup_xrange2,
        float_T = float_T,
        str_in = str_in,
        str_log = str_log
    )
    
    float_deltag = float_g2-float_g1
    float_pka = plm.deltag_to_pka(
        float_deltag = float_deltag,
        float_T = float_T,
        str_in = str_in,
        str_log = str_log
    )

    return float_deltag, float_pka

def get_pka_time(
    str_file: str,
    tup_xrange1: tuple,
    tup_xrange2: tuple,
    str_save: str = None,
    float_T: float = None,
    str_in: str = '../plm.in',
    str_log: str = '../plm.log'
):

    list_file = glob.glob( f'*{str_file}' )
    int_nfile = len(list_file)
    np_pka = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('pka', 'f4')])
    np_deltag = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('deltag', 'f4')])
    for int_i in range(int_nfile-1):
        float_deltag, float_pka = get_pka(
            str_file = f'analysis.{int_i}.{str_file}',
            tup_xrange1 = tup_xrange1,
            tup_xrange2 = tup_xrange2,
            float_T = float_T,
            str_in = str_in,
            str_log = str_log
        )
        np_pka[int_i] = (int_i+1, float_pka)
        np_deltag[int_i] = (int_i+1, float_deltag)
    float_deltag, float_pka = get_pka(
        str_file = str_file,
        tup_xrange1 = tup_xrange1,
        tup_xrange2 = tup_xrange2,
        float_T = float_T,
        str_in = str_in,
        str_log = str_log
    )
    np_pka[int_nfile-1] = (int_nfile, float_pka)
    np_deltag[int_nfile-1] = (int_nfile, float_deltag)
    if str_save:
        np.savetxt(f'{str_save}_pka.csv', np_pka, header=' '.join(np_pka.dtype.names))
        np.savetxt(f'{str_save}_deltag.csv', np_deltag, header=' '.join(np_deltag.dtype.names))

get_pka_time(
    str_file = 'dist_vp_c_fes.grid',
    tup_xrange1 = (1.16,1.66),
    tup_xrange2 = (2.56,'1M'),
    str_save = 'dist_vp_c'
)

'''
run(
    str_file = 'dist_vp_c_fes.grid',
    tup_xrange1 = (1.16,1.66),
    tup_xrange2 = (2.56,'1M')
)

run(
    str_file = 'dist_vp_o_1_2_fes.grid',
    tup_xrange1 = (1.0,1.3),
    tup_xrange2 = (2.4,8.05)
)
'''
