import dpdata
from tf_dpmd_kit import model_devi
import os

dp_sys = dpdata.System(file_name="916100.lammpstrj", fmt='lammps/dump', type_map=["O","H","C"])

str_dir_iter = "../../"

model_devi.def_model_devi_atom(
    list_model = [
        os.path.join(str_dir_iter, "00.train/000/frozen_model.pb"),
        os.path.join(str_dir_iter, "00.train/001/frozen_model.pb"),
        os.path.join(str_dir_iter, "00.train/002/frozen_model.pb"),
        os.path.join(str_dir_iter, "00.train/003/frozen_model.pb"),
    ],
    dp_sys = dp_sys,
    float_lowthred = 0.15,
    #str_save = "devi_f_atom.pdf",
    #float_ylim = 2.6,
    )
