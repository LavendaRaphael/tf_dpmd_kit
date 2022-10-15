import subprocess
import os

int_temperature = 330

def T2KbT(
    float_T: float,
) -> float:

    float_Avogadro = 6.02214076e23
    # J*K^-1
    float_Kb = 1.380649e-23
    # KJ*mol^-1
    float_KbT = float_Kb*float_Avogadro*float_T/1000.0
    return float_KbT

float_KbT = T2KbT(int_temperature)

str_cmd = f'source {os.environ["HOME"]}/.config/.tianff;'
str_cmd += 'source ${homedir}/.local/bin/bashrc_plm.sh;'
str_cmd += f'plumed sum_hills --hills ../../HILLS --stride 20000 --outfile fes.dist_vp_c. --min 0 --max 14 --bin 400 --negbias'

subprocess_results = subprocess.run( str_cmd, shell=True, check=True, text=True, executable='/bin/bash')

print(subprocess_results.stdout)
