import MDAnalysis as mda
from tf_dpmd_kit import analysis
import psutil

dict_sel = {
    'o_w': 'index 0:125',
    'o_0_1_2': 'index 126 127 128',
    'h': 'index 129:382',
    'c': 'index 383'
}

def run():
    mda_u = mda.Universe('traj.lammpstrj', topology_format="LAMMPSDUMP", format="LAMMPSDUMP")
    mda_u.select_atoms("type 1").types = 'O'
    mda_u.select_atoms("type 2").types = 'H'
    mda_u.select_atoms("type 3").types = 'C'
    print(mda_u.trajectory)

    cal = analysis.CarbonicCutoff(
        carbonic_c = mda_u.select_atoms( dict_sel['c'] ),
        carbonic_o = mda_u.select_atoms( dict_sel['o_0_1_2'] ),
        atomg_h = mda_u.select_atoms( dict_sel['h'] ),
        water_o = mda_u.select_atoms( dict_sel['o_w'] ),
    )
    cal.run(
        verbose = True,
    )
    str_save= f'carbonic.csv'
    print(str_save)
    cal.df.to_csv(
        str_save,
        index = False,
    )

p=psutil.Process(3890)
p.wait()
run()
