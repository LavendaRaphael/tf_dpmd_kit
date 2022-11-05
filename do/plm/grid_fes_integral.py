import numpy as np
import math
from tf_dpmd_kit import plm
import glob

def fes_integral(
    str_file: str,
    tup_xrange: tuple,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
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

def fes_1M_correct(
    float_volume: float,
    float_T: float
) -> float:

    float_Avogadro = 6.02214076e23
    float_volume_1M = 1e27/float_Avogadro
    float_KbT = plm.T2KbT(float_T)
    float_correct = -float_KbT*math.log(float_volume_1M/float_volume)

    return float_correct

def get_pka(
    str_file: str,
    tup_xrange1: tuple,
    tup_xrange2: tuple,
    float_volume: float,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
) -> (float, float):

    if not float_T:
        float_T, float_KbT = plm.get_temperature( str_in, str_log )

    float_g1 = fes_integral(
        str_file = str_file,
        tup_xrange = tup_xrange1,
        float_T = float_T,
    )
    
    float_g2 = fes_integral(
        str_file = str_file,
        tup_xrange = tup_xrange2,
        float_T = float_T,
    )

    float_correct = fes_1M_correct(
        float_volume = float_volume,
        float_T = float_T
    )

    float_deltag = float_g2 - float_g1 + float_correct
    float_pka = plm.deltag_to_pka(
        float_deltag = float_deltag,
        float_T = float_T,
    )

    return float_deltag, float_pka

def get_pka_time(
    dict_file: dict,
    tup_xrange1: tuple,
    tup_xrange2: tuple,
    float_volume: float,
    str_save: str = None,
    float_T: float = None,
    str_in: str = None,
    str_log: str = None
):

    int_nfile = len(dict_file)
    np_pka = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('pka', 'f4')])
    np_deltag = np.zeros(shape=(int_nfile), dtype=[('time', 'f4'), ('deltag', 'f4')])
    if not float_T:
        float_T, float_KbT = plm.get_temperature( str_in, str_log )
    for int_i,int_key in enumerate(dict_file):
        float_deltag, float_pka = get_pka(
            str_file = dict_file[int_key],
            tup_xrange1 = tup_xrange1,
            tup_xrange2 = tup_xrange2,
            float_T = float_T,
            float_volume = float_volume
        )
        np_pka[int_i] = (int_key, float_pka)
        np_deltag[int_i] = (int_key, float_deltag)

    if str_save:
        str_save_pka = f'{str_save}_pka.csv'
        print(str_save_pka)
        np.savetxt(str_save_pka, np_pka, header=' '.join(np_pka.dtype.names))
        str_save_deltag = f'{str_save}_deltag.csv'
        print(str_save_deltag)
        np.savetxt(str_save_deltag, np_deltag, header=' '.join(np_deltag.dtype.names))

str_tmp = 'dist_vp_o_1_2_fes'
get_pka_time(
    dict_file = {
        1: f'{str_tmp}.0.dat',
        2: f'{str_tmp}.1.dat',
        3: f'{str_tmp}.2.dat',
        4: f'{str_tmp}.3.dat',
        5: f'{str_tmp}.4.dat',
        6: f'{str_tmp}.5.dat',
        7: f'{str_tmp}.6.dat',
        8: f'{str_tmp}.7.dat',
        9: f'{str_tmp}.8.dat',
        10: f'{str_tmp}.9.dat',
    },
    tup_xrange1 = (0,2),
    tup_xrange2 = (2,14),
    str_save = str_tmp,
    str_in = 'plm.in',
    float_volume = 15.6793091675**3
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
