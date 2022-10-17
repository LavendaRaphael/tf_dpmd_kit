import subprocess
import os

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

    return float_T2KbT

def run(
    str_cv: str,
    str_in: str = '../../plm.in',
    str_log: str = '../../plm.log',
    str_hills: str = '../../HILLS'
):

    float_KbT = get_temperature(str_in, str_log)

    str_cmd = f'source {os.environ["HOME"]}/.config/.tianff &&'
    str_cmd += 'source ${homedir}/.local/bin/bashrc_plm.sh ;'
    str_cmd += f'plumed sum_hills --hills {str_hills} --stride 20000 --outfile fes.{str_cv}. --min 0,-pi --max 14,pi --bin 400,400 --negbias --idw {str_cv} --kt {float_KbT}'
    
    subprocess_results = subprocess.run( str_cmd, shell=True, check=True, text=True, executable='/bin/bash')

    print(subprocess_results.stdout)

run('dist_vp_o_1_2')
