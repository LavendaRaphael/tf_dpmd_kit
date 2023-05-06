from tf_dpmd_kit import convert
import os
import MDAnalysis as mda

def main():

    homedir = os.environ['homedir']

    list_snap = [1]

    mda_u = mda.Universe(
        'lmp.data',
        [
            'test.lammpstrj',
        ],
        topology_format="DATA", format="LAMMPSDUMP"
    )

    convert.mda_to_snap(
        mda_u,
        list_snap
    )

main()

