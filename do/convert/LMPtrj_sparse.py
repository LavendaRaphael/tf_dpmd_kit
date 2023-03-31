from tf_dpmd_kit import convert

def main():

    convert.lmptrj_sparse(
        file_old = 'traj.lammpstrj',
        file_new = 'traj.lammpstrj.new',
        interval = 20,
    )

main()
