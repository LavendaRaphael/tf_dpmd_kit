import os
import numpy as np

def pka_collect(
    list_temperature: list,
    str_data_deltag: str,
    str_data_pka: str,
    str_save_deltag: str,
    str_save_pka: str,
    str_bk: str = ''
) -> None:
    
    int_ntemp = len(list_temp)
    np_pka = np.zeros(int_ntemp, dtype=[('Temperature(K)', 'i4'), ('pKa', 'f4')])
    np_deltag = np.zeros(int_ntemp, dtype=[('Temperature(K)', 'i4'), ('DeltaG(kJ/mol)', 'f4')])
    
    for int_i in range(int_ntemp):
        int_temp = list_temp[int_i]
        str_dir = os.path.join( f'{int_temp}K{str_bk}', 'reweight' )
        with open( os.path.join(str_dir, str_data_pka), 'r' ) as fp:
            float_pka = float(fp.readlines()[-1].split()[1])
        np_pka[int_i] = (int_temp, float_pka)
        with open( os.path.join(str_dir, str_data_deltag), 'r' ) as fp:
            float_deltag = float(fp.readlines()[-1].split()[1])
        np_deltag[int_i] = (int_temp, float_deltag)
    
    np.savetxt(str_save_pka, np_pka, header=' '.join(np_pka.dtype.names))
    np.savetxt(str_save_deltag, np_deltag, header=' '.join(np_deltag.dtype.names))

list_temp = [280,290,300,310,320]

str_tmp = 'dist_vp_o_1_2_fes'

pka_collect(
    list_temperature = list_temp,
    str_data_deltag = f'{str_tmp}.deltag.csv',
    str_data_pka = f'{str_tmp}.pka.csv',
    str_save_deltag = f'reweight.{str_tmp}.deltag.csv',
    str_save_pka = f'reweight.{str_tmp}.pka.csv'
)
