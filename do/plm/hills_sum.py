import subprocess
import os
from tf_dpmd_kit import plm

def run(
    str_cv: str,
    str_in: str = 'plm.in',
    str_log: str = 'plm.log',
    str_hills: str = 'HILLS'
):

    float_T, float_KbT = plm.get_temperature(str_in, str_log)

    str_cmd = f'source {os.environ["HOME"]}/.config/.tianff &&'
    str_cmd += 'source ${homedir}/.local/bin/bashrc_plm.sh ;'
    str_cmd += f'plumed sum_hills --hills {str_hills} --stride 20000 --outfile {str_cv}_fes. --min 0,-pi --max 14,pi --bin 400,400 --negbias --idw {str_cv} --kt {float_KbT}'
    
    subprocess_results = subprocess.run( str_cmd, shell=True, check=True, text=True, executable='/bin/bash')

    print(subprocess_results.stdout)

run(
    str_cv = 'dist_vp_o_1_2',
    str_log = None
)
