import ase.io
import numpy as np
from tf_dpmd_kit import convert

list_snap = [0]

aselist = convert.lmptrj_to_ase(
    list_traj=[
        'traj.lammpstrj'
    ]
)
aselist = convert.aselist_select(
    aselist,
    list_snap
)
conver.aselist_snap(
    aselist,
    list_snap
)
