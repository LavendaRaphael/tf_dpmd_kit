import MDAnalysis as mda
import tf_dpmd_kit.analysis

float_dt = 0.005

mda_u = mda.Universe('traj.lammpstrj', format="LAMMPSDUMP", dt=float_dt)
print(mda_u.trajectory)

tf_dpmd_kit.analysis.msd(
    mda_u = mda_u,
    str_select = 'type 1',
    float_dt = float_dt,
    str_save = 'msd.csv',
    tup_snaprange = (0,15000),
)
