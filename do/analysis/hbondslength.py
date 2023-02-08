import MDAnalysis as mda
from tf_dpmd_kit import analysis

dict_sel = {
    'o_w': 'index 0:125',
    'o_0_2': 'index 126 128',
    'o_0': 'index 126',
    'o_2': 'index 128',
    'o_1': 'index 127',
    'h_w': 'index 129:380',
    'h_0_1': 'index 381 382',
    'h_0': 'index 381',
    'h_1': 'index 382',
    'c': 'index 383'
}

def run(
    dict_pair,
    dict_sel,
):
    mda_u = mda.Universe('../../lmp.data', ['traj.lammpstrj.0','traj.lammpstrj.1'], topology_format="DATA", format="LAMMPSDUMP", dt=0.005)
    mda_u.select_atoms("type 1").types = 'O'
    mda_u.select_atoms("type 2").types = 'H'
    mda_u.select_atoms("type 3").types = 'C'
    print(mda_u.trajectory)

    for str_key, tup_pair in dict_pair.items():
        cal = analysis.HbondLength(
            hydroxyl_o = mda_u.select_atoms( dict_sel[tup_pair[0]] ),
            hydroxyl_h = mda_u.select_atoms( dict_sel[tup_pair[1]] ),
            water_o = mda_u.select_atoms( dict_sel[tup_pair[2]] ),
        )
        cal.run(
            verbose = True,
        )
        str_save= f'hbondslength.{str_key}.csv'
        print(str_save)
        cal.df.to_csv(
            str_save,
            sep = ' ',
            index = False,
        )

dict_pair = {
    'o_0.o_w': ('o_0', 'h_0', 'o_w'),
    'o_2.o_w': ('o_2', 'h_1', 'o_w')
}

run(
    dict_pair = dict_pair,
    dict_sel = dict_sel,
)
