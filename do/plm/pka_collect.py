import os
import numpy as np

list_temp = [290, 310, 330, 350, 370, 390]

str_data_deltag = 'dist_vp_c_deltag.csv'
str_data_pka = 'dist_vp_c_pka.csv'

str_save_deltag = 'reweight_bias.dist_vp_c_deltag.csv'
str_save_pka = 'reweight_bias.dist_vp_c_pka.csv'

int_ntemp = len(list_temp)
np_pka = np.zeros(int_ntemp, dtype=[('Temperature(K)', 'i4'), ('pka', 'f4')])
np_deltag = np.zeros(int_ntemp, dtype=[('Temperature(K)', 'i4'), ('DeltaG(kJ/mol)', 'f4')])

for int_i in range(int_ntemp):
    int_temp = list_temp[int_i]
    if not os.path.exists( f'{int_temp}K' ):
        continue
    with open( os.path.join(f'{int_temp}K', 'reweight_bias', str_data_pka), 'r' ) as fp:
        float_pka = float(fp.readlines()[-1].split()[1])
    np_pka[int_i] = (int_temp, float_pka)
    with open( os.path.join(f'{int_temp}K', 'reweight_bias', str_data_deltag), 'r' ) as fp:
        float_deltag = float(fp.readlines()[-1].split()[1])
    np_deltag[int_i] = (int_temp, float_deltag)

np.savetxt(str_save_pka, np_pka, header=' '.join(np_pka.dtype.names))
np.savetxt(str_save_deltag, np_deltag, header=' '.join(np_deltag.dtype.names))
