import math

def T2KbT(
    float_T: float,
) -> float:

    float_Avogadro = 6.02214076e23
    # J*K^-1
    float_Kb = 1.380649e-23
    # KJ*mol^-1
    float_KbT = float_Kb*float_Avogadro*float_T/1000.0
    return float_KbT

def get_temperature(
    str_in: str = '../../plm.in',
    str_log: str = '../../plm.log'
) -> float:
    
    with open(str_in, 'r') as fp:
        for str_line in fp:
            if 'TEMP=' in str_line:
                list_line = str_line.split()
                break
        for str_key in list_line:
            if 'TEMP=' == str_key[:5]:
                float_T = int(str_key[5:])
                break
    print(float_T)
    float_T2KbT = T2KbT(float_T)

    with open(str_log, 'r') as fp:
        for str_line in fp:
            if 'KbT' in str_line:
                list_line = str_line.split()
                break
        float_KbT = float(list_line[2])
    print(float_KbT)

    if float_T2KbT-float_KbT > 1e4:
        raise

    return float_T, float_T2KbT

def deltag_to_pka(
    float_deltag: float,
    float_T: float = None,
    str_plmin: str = '../../plm.in',
    str_plmlog: str = '../../plm.log'
) -> float:

    if float_T is None:
        float_T, float_KbT = get_temperature(str_plmin, str_plmlog)
    else:
        float_KbT = T2KbT(float_T)

    return float_deltag/(float_KbT*math.log(10))

def pka_to_deltag(
    float_pka: float,
    float_T: float = None,
    str_plmin: str = '../../plm.in',
    str_plmlog: str = '../../plm.log'
) -> float:

    if float_T is None:
        float_T, float_KbT = get_temperature(str_plmin, str_plmlog)
    else:
        float_KbT = T2KbT(float_T)

    return float_pka*(float_KbT*math.log(10))

