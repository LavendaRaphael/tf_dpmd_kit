import pandas as pd
import numpy as np
import time

def run():
    start = time.time()

    df_data = pd.read_csv('carbonic.csv')
    print(df_data)

    df_new = df_data.apply(lambda x: state(x['ncarbonyl'], x['nhydroxyl'], x['dihedral0(rad)'], x['dihedral1(rad)']), axis=1, result_type='expand')
    df_new.columns = ['CO3','0.5','HCO3','1.5','TT','CT','CC','2.5','H3CO3']
    df_new.insert(0, 'frame', df_data['frame'])
    print(df_new)

    df_new.to_csv('carbonic_state.csv', index=False)

    end = time.time()
    print(end-start)

def state(
    ncarbonyl, 
    nhydroxyl,
    alpha,
    beta,
):
    if ncarbonyl == 3:
        # 300 CO3
        int_state = 0
    elif ncarbonyl == 2:
        if nhydroxyl != 1:
            # 201
            int_state = 1
        else:
            # 210 HCO3
            int_state = 2
    elif ncarbonyl == 1:
        if nhydroxyl != 2:
            # 102 111
            int_state = 3
        else:
            # 120 H2CO3
            int_state = conformer(alpha, beta)
    elif ncarbonyl == 0:
        if nhydroxyl != 3:
            # 003 012 021
            int_state = 7
        else:
            # 003 H3CO3
            int_state = 8
    list_re = [None]*9
    list_re[int_state] = 1
    return list_re

def conformer(
    alpha,
    beta
):
    pio2 = np.pi/2
    bool_alpha = ((alpha > -pio2) & (alpha < pio2))
    bool_beta = ((beta > -pio2) & (beta < pio2))

    if bool_alpha:
        if bool_beta:
            # TT
            int_state = 4
        else:
            # TC
            int_state = 5
    else:
        if bool_beta:
            # CT
            int_state = 5
        else:
            # CC
            int_state = 6

    return int_state

run()
